import sys
sys.path.insert(0, "../")
import networkx as nx
from networkx.algorithms.components import connected_components
from networkx.algorithms.community.centrality import girvan_newman
from networkx.algorithms.community.quality import modularity
import os
from src.core.network_builder import build_network

def create_slices(network_file, output_file_name):
    G = build_network(network_file)
    G.remove_nodes_from(list(nx.isolates(G)))
    print("after removing orphans - number of edges: {}, nodes: {}".format(len(G.edges),
                                                                                    len(G.nodes)))

    reached_optimal_modularity = False
    optimized_connected_components = girvan_newman(G)
    optimal_components=[]
    n_modules=[]
    n_large_modules=[]
    modularity_scores=[]
    optimal_modularity=-1
    while not reached_optimal_modularity:

        try:
            cur_components = sorted(next(optimized_connected_components))
        except StopIteration:
            break

        cur_modularity = modularity(G, cur_components, weight='weight')
        if cur_modularity > optimal_modularity or True:
            optimal_modularity = cur_modularity
            optimal_components = cur_components

            edges_to_remove = []
            for cur_edge in G.edges:
                included = False
                for cur_cc in optimal_components:
                    if cur_edge[0] in cur_cc and cur_edge[1] in cur_cc:
                        included = True
                if not included:
                    edges_to_remove.append(cur_edge)

            G.remove_edges_from(edges_to_remove)

            n_modules.append(len(cur_components))
            n_large_modules.append(len([a for a in cur_components if len(a) > 3]))
            modularity_scores.append(optimal_modularity)
            print("new: ", len(cur_components), len([a for a in cur_components if len(a) > 3]), optimal_modularity)


        with open(output_file_name, 'w+') as f:
            f.write(str(n_modules)+"\n")
            f.write(str(n_large_modules)+"\n")
            f.write(str([a*100 for a in modularity_scores])+"\n")

    print("modularity: ", modularity(G, list([G.subgraph(c) for c in connected_components(G)]), weight='weight'))


def read_preprocessed_slices(file_path):
    modules = []

    with open(file_path, 'r') as f:
        line = f.readline()
        while line != "":
            line = f.readline()
            if line.startswith("cc"):
                modules.append(f.readline().strip()[1:-1].split(', '))

        return modules
