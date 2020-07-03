import networkx as nx
import pandas as pd


def build_network(network_file):
    """"""
    edges_dataset = pd.read_csv(network_file, sep='\t', header=None, dtype=str)
    edges = []
    for ind, row in edges_dataset.iterrows():
        # if row.iloc[0]!=row.iloc[2]:
        edges.append((row.iloc[0], row.iloc[2]))
    G = nx.Graph()
    G.add_edges_from(edges)
    nx.set_node_attributes(G, 0, 'score')

    return G
