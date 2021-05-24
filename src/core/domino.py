import sys
sys.path.insert(0, "../")

import random
import os

import pandas as pd
import numpy as np
import pickle
import multiprocessing
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

import pcst_fast
import networkx as nx
from networkx.algorithms.community.quality import modularity
from networkx.algorithms.community.centrality import girvan_newman
from networkx.algorithms.components import connected_components

from functools import reduce
from src.utils.graph_influence_linear_th import linear_threshold
from src.core.preprocess_slices import read_preprocessed_slices
from src.core.network_builder import build_network
import src.constants as constants

G_modularity = None


def extract_scores(scores_file):
    """"""
    scores = pd.read_csv(scores_file, sep='\t', index_col=0, header=None, dtype=str)
    if "pval" in scores.columns:
        scores["score"] = scores["pval"]
    else:
        scores["score"] = 1
    return scores


def add_scores_to_nodes(G, scores):
    """"""
    inds = []
    for nd in G.nodes:
        G.nodes[nd]["pertubed_node"] = False
        G.nodes[nd]["score"] = 0

    for ind, row in scores.iterrows():
        if ind in G.nodes:
            inds.append(ind)
            G.nodes[ind]["score"] = row["score"]
            G.nodes[ind]["pertubed_node"] = row["score"] > 0  # binarizing the activeness

    return G


def create_subgraph(params):
    cur_module = params
    global G_modularity
    nodes = set(cur_module)
    res = G_modularity.subgraph(list(nodes))
    return res


def prune_network_by_modularity(G, modules, cache_file):
    global G_modularity
    if os.path.exists(cache_file) and constants.USE_CACHE:
        print(f'fetch cache file for subnetworks {cache_file}')
        G_modularity = pickle.load(open(cache_file, 'rb'))
        for n in G_modularity:
            G_modularity.nodes[n]['pertubed_node'] = G.nodes[n]['pertubed_node']
        print('pkl is loaded')
        return G_modularity

    print(f'generating subgraphs...')
    G_modularity = G
    print(
        f"Before slicing: n of cc:{len(list(connected_components(G_modularity)))}, n of nodes: {len(G_modularity.nodes)}, n of edges, {len(G_modularity.edges)}")
    p = multiprocessing.Pool(constants.N_OF_THREADS)

    G_modules = p.map(create_subgraph, [m for m in modules])
    p.close()
    # print(f'{modules}')
    print(f'# of modules after extraction: {len(G_modules)}')
    G_modularity = nx.algorithms.operators.union_all(G_modules)
    print(
        f"After slicing: n of cc:{len(list(connected_components(G_modularity)))}, n of nodes: {len(G_modularity.nodes)}, n of edges, {len(G_modularity.edges)}")
    pickle.dump(G_modularity, open(cache_file, 'wb+'))
    print('subgraphs\' pkl is saved')


def prune_network_by_modularity_old(G, modules, dummy):
    G_modularity = G.copy()
    edges_to_remove = []
    for cur_edge in G_modularity.copy().edges:
        in_cc = False
        for cur_module in modules:
            if cur_edge[0] in cur_module and cur_edge[1] in cur_module:
                in_cc = True
        if not in_cc:
            edges_to_remove.append(cur_edge)

    G_modularity.remove_edges_from(edges_to_remove)
    return G_modularity


def get_pcst_prize(G_cc, prize_factor, n_steps):
    prizes = {}
    p_cc = linear_threshold(G_cc, [n for n in G_cc.nodes if G_cc.nodes[n]['pertubed_node'] > 0], steps=n_steps)
    for p_node in G_cc.nodes:
        prizes[p_node] = 0
    for i_cur_layer, cur_layer in enumerate(p_cc):
        for cur_node in cur_layer:
            prizes[cur_node] += prize_factor ** i_cur_layer

    return prizes

def run_pcst(G_cc, i_cc, labels, n_steps, nodes, prize_factor):
    ## set prize ##
    prizes = get_pcst_prize(G_cc, prize_factor, n_steps)
    vertices_prizes = []
    for cur_node in nodes:
        vertices_prizes.append(
            G_cc.nodes[cur_node]["pertubed_node"] if G_cc.nodes[cur_node]["pertubed_node"] else prizes[cur_node])

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
    cur_components = [G_optimized.subgraph(c) for c in connected_components(G_optimized)]
    cur_modularity = modularity(G_optimized, cur_components, weight='weight')
    if cur_modularity >= modularity_score_objective:
        return True, best_modularity

        if len(n_nodes) < 4:
            G_optimized.remove_nodes_from(n_nodes)

    cur_components = [G_optimized.subgraph(c) for c in connected_components(G_optimized)]
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


def get_putative_modules(G, full_G=None, improvement_delta=0, modularity_score_objective=1, module_threshold=0.05,
                         n_cc=1.0):
    """"""

    if full_G == None:
        full_G = G
    G_optimized = G.copy()

    # clean subslice from cycles and isolated nodes
    G_optimized.remove_edges_from(list(nx.selfloop_edges(G_optimized)))
    G_optimized.remove_nodes_from(list(nx.isolates(G_optimized)))

    # check subslice enrichment for active nodes
    pertubed_nodes = [cur_node for cur_node in full_G.nodes if full_G.nodes[cur_node]["pertubed_node"]]
    pertubed_nodes_in_cc = [n for n in G_optimized.nodes if G_optimized.nodes[n]["pertubed_node"]]
    n_nodes = list(G_optimized.nodes)
    sig_score = hypergeom.sf(len(pertubed_nodes_in_cc), len(full_G.nodes), len(pertubed_nodes),
                             len(n_nodes)) \
                + hypergeom.pmf(len(pertubed_nodes_in_cc), len(full_G.nodes), len(pertubed_nodes),
                                len(n_nodes))

    sig_score = sig_score / n_cc

    # if subslice is not enriched for active nodes split in into putative modules. otherwise, report it as a single putative module
    # print(f'{sig_score}<{module_threshold} and {len(G_optimized.nodes)}<30')
    is_enriched_sublice = (len(G_optimized.nodes) < 100) or len(
        G_optimized.nodes) == 0  # sig_score<module_threshold and l

    break_loop = is_enriched_sublice
    best_modularity = -1
    while not break_loop:
        break_loop, best_modularity = split_subslice_into_putative_modules(G_optimized, improvement_delta,
                                                                           modularity_score_objective, best_modularity)

    G_optimized.remove_nodes_from(list(nx.isolates(G_optimized)))

    cc_optimized = [] if len(G_optimized.nodes) == 0 else [G_optimized.subgraph(c) for c in
                                                           connected_components(G_optimized)]

    return G_optimized, cc_optimized


def retain_relevant_slices(G_original, module_sig_th):
    global G_modularity

    pertubed_nodes = []
    for cur_node in G_modularity.nodes():
        if G_modularity.nodes[cur_node]["pertubed_node"]:
            pertubed_nodes.append(cur_node)

    ccs = [G_modularity.subgraph(c) for c in connected_components(G_modularity)]
    params = []
    p = multiprocessing.Pool(constants.N_OF_THREADS)
    n_G_original = len(G_original)
    n_pertubed_nodes = len(pertubed_nodes)
    pertubed_nodes_in_ccs = []
    print(f"number of slices: {len(list(ccs))}")
    for i_cur_cc, cur_cc in enumerate(ccs):
        pertubed_nodes_in_ccs.append(
            len([cur_node for cur_node in cur_cc if G_modularity.nodes[cur_node]["pertubed_node"]])) 
    perturbation_factor = min(0.7, (float(n_pertubed_nodes) / n_G_original) * (
                1 + 100 / n_G_original ** 0.5))

    for i_cur_cc, cur_cc in enumerate(ccs):
        params.append([n_G_original, cur_cc, i_cur_cc, n_pertubed_nodes, perturbation_factor])

    res = [a for a in p.map(pf_filter, params) if a is not None]
    print(f'# of slices after perturbation TH: {len(res)}/{len(params)}')
    p.close()
    if len(res) == 0:
        return nx.Graph(), [], []
    large_modules, sig_scores = zip(*res)
    fdr_bh_results = fdrcorrection0(sig_scores, alpha=module_sig_th, method='indep',
                                    is_sorted=False)

    # print(fdr_bh_results)
    # print(f'min: {min(list(fdr_bh_results[1]))}')
    passed_modules = [cur_cc for cur_cc, is_passed_th in zip(large_modules, fdr_bh_results[0]) if is_passed_th]
    return nx.algorithms.operators.union_all(passed_modules) if len(passed_modules) > 0 else nx.Graph(), [list(m.nodes)
                                                                                                          for m in
                                                                                                          passed_modules], \
           fdr_bh_results[1]


def pf_filter(params):
    global G_modularity
    n_G_original, cur_cc, i_cur_cc, n_pertubed_nodes, perturbation_factor = params
    pertubed_nodes_in_cc = [cur_node for cur_node in cur_cc if G_modularity.nodes[cur_node]["pertubed_node"]]
    if len(cur_cc) < 4 or n_pertubed_nodes == 0 or not (
            len(pertubed_nodes_in_cc) / float(len(cur_cc)) >= perturbation_factor or len(pertubed_nodes_in_cc) / float(
            n_pertubed_nodes) >= 0.1):
        return None
    else:
        score = hypergeom.sf(len(pertubed_nodes_in_cc), n_G_original, n_pertubed_nodes,
                             len(cur_cc)) \
                + hypergeom.pmf(len(pertubed_nodes_in_cc), n_G_original, n_pertubed_nodes,
                                len(cur_cc))
        return (cur_cc, score)


def analyze_slice(params):
    G, cc, i_cc, n_steps, relevant_slices, prize_factor, module_threshold = params
    G_cc = nx.subgraph(G, cc)
    nodes = list(G_cc.nodes)
    labels = {n: G_cc.nodes[n] for n in nodes}
    n_pertubed_nodes = sum([G.nodes[a]["pertubed_node"] for a in G.nodes])
    prize_factor = max(0, 1 - 3 * n_pertubed_nodes / float(len(G.nodes)))
    # print(f'active gene ratio: {n_pertubed_nodes}/{len(G_cc.nodes)}')
    # print(f"prize factor: {prize_factor}")
    edges, edges_grid = run_pcst(G_cc, i_cc, labels, n_steps, nodes, prize_factor)
    G_subslice = nx.Graph()
    G_subslice.add_edges_from([(nodes[edges_grid[e][0]], nodes[edges_grid[e][1]]) for e in edges])
    nx.set_node_attributes(G_subslice, {n: labels[n] for n in G_subslice.nodes})
    modularity_score_objective = np.log(len(G_subslice.nodes)) / np.log(len(G.nodes)) if len(
        G_subslice.nodes) > 10 else -1
    subslice_after_ng, putative_modules_of_slice = get_putative_modules(G_subslice, G, improvement_delta=10 ** -2,
                                                                        modularity_score_objective=modularity_score_objective,
                                                                        n_cc=len(relevant_slices),
                                                                        module_threshold=module_threshold)

    return putative_modules_of_slice


def get_final_modules(G, G_putative_modules):
    module_sigs = []
    for i_cur_module, cur_G_module in enumerate(G_putative_modules):
        pertubed_nodes_in_cc = [cur_node for cur_node in cur_G_module.nodes if G.nodes[cur_node]["pertubed_node"]]
        pertubed_nodes = [cur_node for cur_node in G.nodes if G.nodes[cur_node]["pertubed_node"]]

        sig_score = hypergeom.sf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes),
                                 len(cur_G_module.nodes)) \
                    + hypergeom.pmf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes),
                                    len(cur_G_module.nodes))

        final_module_threshold = 0.05 / len(G_putative_modules)
        if sig_score <= final_module_threshold:
            module_sigs.append((cur_G_module, sig_score / len(G_putative_modules)))

    module_sigs = sorted(module_sigs, key=lambda a: a[1])
    return [a[0] for a in module_sigs]


def main(active_genes_file, network_file, slices_file=None, slice_threshold=0.3, module_threshold=0.05, prize_factor=0,
         n_steps=20):
    print("start running DOMINO...")
    if os.path.exists(f'{network_file}.pkl') and constants.USE_CACHE:
        G = pickle.load(open(f'{network_file}.pkl', 'rb'))
        print(f'network\' pkl is loaded: {network_file}.pkl')
    else:
        print(f'generating graph from {network_file}')
        G = build_network(network_file)
        pickle.dump(G, open(f'{network_file}.pkl', 'wb+'))
        print(f'network\' pkl is saved: {network_file}.pkl')

    print("done building network")
    # assign activeness to nodes
    scores = extract_scores(active_genes_file)
    G = add_scores_to_nodes(G, scores)

    modularity_connected_components = read_preprocessed_slices(slices_file)

    global G_modularity
    prune_network_by_modularity(G, modularity_connected_components, os.path.join(os.path.split(slices_file)[0],
                                                                                 os.path.split(network_file)[1].split(
                                                                                     ".")[0] + "." +
                                                                                 os.path.split(slices_file)[1].split(
                                                                                     ".")[0] + ".pkl"))
    G_modularity, relevant_slices, qvals = retain_relevant_slices(G, slice_threshold)
    print(f'{len(relevant_slices)} relevant slices were retained with threshold {slice_threshold}')
    params = []
    for i_cc, cc in enumerate(relevant_slices):
        params.append([G, cc, i_cc, n_steps, relevant_slices, prize_factor, module_threshold])
    p = multiprocessing.Pool(constants.N_OF_THREADS)
    putative_modules = reduce(lambda a, b: a + b, p.map(analyze_slice, params), [])
    p.close()
    print(f'n of putative modules: {len(putative_modules)}')
    final_modules = get_final_modules(G, putative_modules)
    print(
        f'n of final modules: {len(final_modules)} (n={[len(list(m)) for m in final_modules]})')
    return final_modules
