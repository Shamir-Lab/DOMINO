import sys
sys.path.insert(0, "../")

import random
import os

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

import pcst_fast
import networkx as nx
from networkx.algorithms.community.quality import modularity
from networkx.algorithms.community.centrality import girvan_newman
from networkx.algorithms.components import connected_components

from src.utils.graph_influence_linear_th import linear_threshold
from src.core.preprocess_slices import read_preprocessed_slices
from src.core.network_builder import build_network


def extract_scores(scores_file):
    """"""
    scores = pd.read_csv(scores_file, sep='\t', index_col=0, header=None)
    if "pval" in scores.columns:
        scores["score"] = scores["pval"]
    else:
        scores["score"] = 1
    return scores

def add_scores_to_nodes(G, scores):
    """"""
    inds = []
    for nd in G.nodes:
        G.nodes[nd]["pertubed_node"]=False
        G.nodes[nd]["score"]=0

    for ind, row in scores.iterrows():
        if ind in G.nodes:
            inds.append(ind)
            G.nodes[ind]["score"] = row["score"]
            G.nodes[ind]["pertubed_node"] = row["score"] > 0  # binarizing the activeness

    return G

def prune_network_by_modularity(G, modules):
    G_modularity=G.copy()
    print(f"Before slicing: n of cc:{len(list(connected_components(G_modularity)))}, n of nodes: {len(G_modularity.nodes)}, n of edges, {len(G_modularity.edges)}")
    edges_to_remove=[]
    for cur_edge in G_modularity.copy().edges:
        in_cc = False
        for cur_module in modules:
            if cur_edge[0] in cur_module and cur_edge[1] in cur_module:
                in_cc=True
        if not in_cc:
            edges_to_remove.append(cur_edge)

    G_modularity.remove_edges_from(edges_to_remove)
    print(f"After slicing: n of cc:{len(list(connected_components(G_modularity)))}, n of nodes: {len(G_modularity.nodes)}, n of edges, {len(G_modularity.edges)}")
    return G_modularity

def get_pcst_prize(G_cc, prize_factor, n_steps):
    prizes={}
    p_cc = linear_threshold(G_cc, [n for n in G_cc.nodes if G_cc.nodes[n]['pertubed_node'] > 0], steps=n_steps)
    for p_node in G_cc.nodes:
        prizes[p_node] = 0
    for i_cur_layer, cur_layer in enumerate(p_cc):
        for cur_node in cur_layer:
            prizes[cur_node] += prize_factor ** i_cur_layer

    return prizes


def get_permuted_scores(G_cc, labels, n_steps, nodes, prize_factor, n_permutations):
    permuted_scores = {}
    for cur_node in nodes:
        permuted_scores[cur_node] = []
    for i_iteration in np.arange(n_permutations):
        c_nodes = list(nodes)
        random.shuffle(c_nodes)
        mapping = {k: v for k, v in zip(nodes, c_nodes)}
        p_cc = nx.relabel_nodes(G_cc, mapping)
        nx.set_node_attributes(p_cc, labels)
        lt_genes = linear_threshold(p_cc, [n for n in p_cc.nodes if p_cc.nodes[n]['pertubed_node'] > 0], steps=n_steps)
        for p_node in G_cc.nodes:
            permuted_scores[p_node] = []
        for i_cur_layer, cur_layer in enumerate(lt_genes):
            for cur_node in cur_layer:
                permuted_scores[cur_node].append(prize_factor ** i_cur_layer)
    return permuted_scores


def run_pcst(G_cc, i_cc, labels, n_steps, nodes, prize_factor, n_permutations=1000):

    ## set prize ##
    prizes = get_pcst_prize(G_cc, prize_factor, n_steps)
    vertices_prizes = []
    for cur_node in nodes:
        vertices_prizes.append(G_cc.nodes[cur_node]["pertubed_node"] if G_cc.nodes[cur_node]["pertubed_node"] else prizes[cur_node])

    ## set cost ##
    edges_grid = []
    for cur_edge in G_cc.edges:
        edges_grid.append([nodes.index(cur_edge[0]), nodes.index(cur_edge[1])])

    edges_costs = []
    for cur_edge in edges_grid:
        u_score = 0 if G_cc.nodes[nodes[cur_edge[0]]]["pertubed_node"] else 0.9999
        v_score = 0 if G_cc.nodes[nodes[cur_edge[1]]]["pertubed_node"] else 0.9999

        edges_costs.append(np.min([u_score, v_score]))

    ## find pcst component by running pcst fast##
    root = -1
    num_clusters = 1
    pruning = 'strong'  # 'none'
    verbosity_level = 0
    vertices, edges = pcst_fast.pcst_fast(edges_grid, vertices_prizes, edges_costs, root, num_clusters, pruning,
                                          verbosity_level)

    return edges, edges_grid


def split_subslice_into_putative_modules(G_optimized, improvement_delta, modularity_score_objective, best_modularity):

    n_nodes=list(G_optimized.nodes)
    cur_components = [G_optimized.subgraph(c) for c in connected_components(G_optimized)]
    cur_modularity = modularity(G_optimized, cur_components, weight='weight')
    # print("cur_modularity: {}".format(cur_modularity))
    if cur_modularity >= modularity_score_objective:
        return True, best_modularity

        if len(n_nodes) < 4:
            G_optimized.remove_nodes_from(n_nodes)
    # print("after optimizing for sig modules - number of edges: {}, nodes: {}".format(len(G_optimized.edges),
    #                                                                                  len(G_optimized.nodes)))
    cur_components = [G_optimized.subgraph(c) for c in connected_components(G_optimized)]
    # print("number of cc: ", len(cur_components))
    if len(cur_components) == 0:
        return True, best_modularity

    optimized_connected_components = girvan_newman(G_optimized)
    cur_components = sorted(next(optimized_connected_components))
    cur_modularity = modularity(G_optimized, cur_components, weight='weight')
    if cur_modularity <= best_modularity + improvement_delta:
        return True, best_modularity

    else:
        optimal_components = cur_components

        edges_to_remove = []
        for cur_edge in G_optimized.edges:
            included = False
            for n_nodes in optimal_components:
                if cur_edge[0] in n_nodes and cur_edge[1] in n_nodes:
                    included = True
            if not included:
                edges_to_remove.append(cur_edge)

        G_optimized.remove_edges_from(edges_to_remove)

        return False, cur_modularity


def get_putative_modules(G, full_G=None, improvement_delta=0, modularity_score_objective=1, module_threshold=0.05, n_cc=1.0):
    """"""

    if full_G==None:
        full_G=G
    G_optimized = G.copy()

    # clean subslice from cycles and isolated nodes
    G_optimized.remove_edges_from(list(nx.selfloop_edges(G_optimized)))
    G_optimized.remove_nodes_from(list(nx.isolates(G_optimized)))

    # check subslice enrichment for active nodes
    pertubed_nodes = [cur_node for cur_node in full_G.nodes if full_G.nodes[cur_node]["pertubed_node"]]
    pertubed_nodes_in_cc=[n for n in G_optimized.nodes if G_optimized.nodes[n]["pertubed_node"]]
    n_nodes=list(G_optimized.nodes)
    sig_score = hypergeom.sf(len(pertubed_nodes_in_cc), len(full_G.nodes), len(pertubed_nodes),
                             len(n_nodes)) \
                + hypergeom.pmf(len(pertubed_nodes_in_cc), len(full_G.nodes), len(pertubed_nodes),
                                len(n_nodes))

    sig_score=sig_score/n_cc

    # if subslice is not enriched for active nodes split in into putative modules. otherwise, report it as a single putative module
    is_enriched_sublice = (sig_score<module_threshold and len(G_optimized.nodes)<10) or len(G_optimized.nodes)==0

    break_loop = is_enriched_sublice
    best_modularity=-1
    while not break_loop:
        break_loop, best_modularity=split_subslice_into_putative_modules(G_optimized, improvement_delta, modularity_score_objective, best_modularity)

    G_optimized.remove_nodes_from(list(nx.isolates(G_optimized)))

    cc_optimized = [] if len(G_optimized.nodes) == 0 else [G_optimized.subgraph(c) for c in connected_components(G_optimized)]


    return G_optimized, cc_optimized


def retain_relevant_slices(G_original, G_modularity, module_sig_th):

    pertubed_nodes = []
    for cur_node in G_modularity.nodes():
        if G_modularity.nodes[cur_node]["pertubed_node"]:
            pertubed_nodes.append(cur_node)

    ccs=[G_modularity.subgraph(c) for c in connected_components(G_modularity)]
    sig_scores=[]
    passed_modules=[]
    large_modules=[]
    G_optimized=G_modularity.copy()
    perturbation_factor = np.log2(len([n for n in G_original.nodes if G_original.nodes[n]['pertubed_node']]))
    print(f'n_of_ccs: {len(ccs)}')
    for cur_cc in ccs:

        pertubed_nodes_in_cc = [cur_node for cur_node in cur_cc if G_modularity.nodes[cur_node]["pertubed_node"]]
        if len(cur_cc) < 4 or not (len(pertubed_nodes_in_cc) > perturbation_factor or len(pertubed_nodes_in_cc)/float(len(cur_cc))>=0.2) :
            G_optimized.remove_nodes_from(cur_cc)
        else:
            large_modules.append(cur_cc)

            sig_scores.append(hypergeom.sf(len(pertubed_nodes_in_cc), len(G_original.nodes), len(pertubed_nodes),
                                           len(cur_cc)) \
                        + hypergeom.pmf(len(pertubed_nodes_in_cc), len(G_original.nodes), len(pertubed_nodes),
                                        len(cur_cc)))
            print(f'HG({len(pertubed_nodes_in_cc)}, {len(G_original.nodes)}, {len(pertubed_nodes)}, {len(cur_cc)}) = {sig_scores[-1]}')

    fdr_bh_results = fdrcorrection0(sig_scores, alpha=module_sig_th, method='indep',
                                    is_sorted=False)

    print(fdr_bh_results[1])
    i=0
    for cur_cc, is_passed_th, qval in zip(large_modules, fdr_bh_results[0], fdr_bh_results[1]):
        if is_passed_th:
            passed_modules.append(list(cur_cc.nodes))
        else:
            G_optimized.remove_nodes_from(cur_cc)

        i += 1

    return G_optimized, passed_modules, fdr_bh_results[1]


def analyze_slice(G, cc, i_cc, n_steps, relevant_slices, prize_factor, module_threshold, n_permutations):
    G_cc = nx.subgraph(G, cc)
    nodes = list(G_cc.nodes)
    labels = {n: G_cc.nodes[n] for n in nodes}
    edges, edges_grid = run_pcst(G_cc, i_cc, labels, n_steps, nodes, prize_factor, n_permutations)
    G_subslice = nx.Graph()
    G_subslice.add_edges_from([(nodes[edges_grid[e][0]], nodes[edges_grid[e][1]]) for e in edges])
    nx.set_node_attributes(G_subslice, {n: labels[n] for n in G_subslice.nodes})
    modularity_score_objective = np.log(len(G_subslice.nodes)) / np.log(len(G.nodes)) if len(G_subslice.nodes) > 10 else -1
    subslice_after_ng, putative_modules_of_slice = get_putative_modules(G_subslice, G, improvement_delta=10 ** -2,
                                                                        modularity_score_objective=modularity_score_objective,
                                                                        n_cc=len(relevant_slices), module_threshold=module_threshold)

    return [list(sg.nodes) for sg in putative_modules_of_slice]

def get_final_modules(G, putative_modules):
    module_sigs = []
    for i_cur_module, cur_module in enumerate(putative_modules):
        pertubed_nodes_in_cc = [cur_node for cur_node in cur_module if G.nodes[cur_node]["pertubed_node"]]
        pertubed_nodes = [cur_node for cur_node in G.nodes if G.nodes[cur_node]["pertubed_node"]]

        sig_score = hypergeom.sf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes),
                                 len(cur_module)) \
                    + hypergeom.pmf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes),
                                    len(cur_module))

        if sig_score <= 0.05 / len(putative_modules):
            module_sigs.append((cur_module, sig_score / len(putative_modules)))


    module_sigs = sorted(module_sigs, key=lambda a: a[1])
    return [a[0] for a in module_sigs]

def main(active_genes_file, network_file, slices_file=None, slice_threshold=0.3, module_threshold=0.05, prize_factor=0.7, n_steps=20, n_permutations=1000):

    network_file_name = os.path.join(network_file)

    G = build_network(network_file_name)

    # assign activeness to nodes
    scores = extract_scores(active_genes_file)
    G = add_scores_to_nodes(G, scores)

    modularity_connected_components = read_preprocessed_slices(slices_file)

    G_modularity=prune_network_by_modularity(G, modularity_connected_components)
    G_modularity, relevant_slices, qvals = retain_relevant_slices(G, G_modularity, slice_threshold)

    putative_modules = []
    for i_cc, cc in enumerate(relevant_slices):
        putative_modules += analyze_slice(G, cc, i_cc, n_steps, relevant_slices, prize_factor, module_threshold, n_permutations)

    final_modules=get_final_modules(G, putative_modules)
    return get_final_modules(G, putative_modules)