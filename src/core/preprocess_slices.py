import sys
sys.path.insert(0, "../")
import networkx as nx
import pandas as pd
import community as community_louvain
import numpy as np

def create_slices(network_file, output_file_name, resolution=0.15):
    df = pd.read_csv(network_file, sep='\t')
    df.columns = ["node_1", "edge_type", "node_2"]
    G = nx.from_pandas_edgelist(df, 'node_1', 'node_2')

    partition = community_louvain.best_partition(G, resolution=resolution, random_state=1)  # 0.1
    prt = {k: [] for k in np.arange(len(np.unique(list(partition.values()))))}
    for k, v in partition.items():
        prt[v].append(k)

    i = 0
    with open(output_file_name, 'w+') as f:
        f.write(f'# of cc after modularity optimization: {len(prt.keys())}\n')
        for k, v in prt.items():
            if len(v) >= 10:
                f.write(f'cc #{i}: n={len(v)}\n[{", ".join(v)}]\n')
                i += 1


def read_preprocessed_slices(file_path):
    modules = []

    with open(file_path, 'r') as f:
        line = f.readline()
        while line != "":
            line = f.readline()
            if line.startswith("cc"):
                modules.append(f.readline().strip()[1:-1].split(', '))

        return modules
