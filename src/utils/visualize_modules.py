
import sys
sys.path.insert(0, '../../')
import numpy as np
import os
import time
import shutil
import json
import pandas as pd

from src import constants
from src.utils.scripts import format_script
from src.utils.ensembl2gene_symbol import  e2g_convertor
import zipfile

import multiprocessing
from functools import reduce
SH_MODULE_NAME = "module"
SH_NUM_GENES = "#_genes"
SH_ENRICHED = "enriched_groups"
SH_DETAILS = "more_details"

SH_TABLE_HEADERS = [SH_MODULE_NAME, SH_NUM_GENES, SH_ENRICHED, SH_DETAILS]

MODULE_TH = 10

def zipdir(path_to_zip, zip_file_path):
    ziph = zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(path_to_zip):
        for file in files:
            ziph.write(os.path.join(root, file))

def get_network_genes(network_file_name):
    network_df = pd.read_csv(network_file_name, sep="\t")
    src = np.array(network_df.loc[:, 0])
    dst = np.array(network_df.loc[:, 2])
    vertices = list(set(np.append(src, dst)))
    return vertices

def remove_subgraph_self_loops(nodes_to_remove, network_file_name):
    if len(nodes_to_remove) == 0:
        return network_file_name
    network_df = pd.read_csv(network_file_name, sep="\t")
    filtered_network = network_df[network_df.loc[:,0]!=network_df.loc[:,2]]
    new_file_name = os.path.splitext(network_file_name) + "_no_loops" +".sif"
    filtered_network.to_csv(new_file_name, sep="\t", index=False)
    return filtered_network

def remove_subgraph_by_nodes(nodes_to_remove, network_file_name, ts=str(time.time())):
    if len(nodes_to_remove) == 0:
        return network_file_name
    network_df = pd.read_csv(network_file_name, sep="\t")
    filtered_network = network_df[~(network_df.loc[:,0].isin(nodes_to_remove) | network_df.loc[:,2].isin(nodes_to_remove))]
    new_file_name = os.path.splitext(network_file_name)[0] + ts +".sif"
    filtered_network.to_csv(new_file_name, sep="\t", index=False)
    return new_file_name


def summary_intergrative_reports(all_hg_reports, modules_summary, total_hg_report, algo_name, module_genes, report_file_name, dataset_name):
    general_algo_report(algo_name, all_hg_reports, module_genes, modules_summary, report_file_name, total_hg_report, dataset_name)


def disease_algo_report(algo_name, disease_name, expected_genes, module_genes, modules_summary, report_file_name, dataset_name):

    disease_data = {
        "disease_name": disease_name,
        "num_of_modules": len(modules_summary),
        "TP+FN_(_true_)": len(expected_genes),
        "TP+TN_(_retrieved_)": len(module_genes),
        "TP/(TP+TN)_(_precision_)": 0,
        "TP/(TP+FN)_(_recall_)": 0,
        "F1": 0,
        "TP": 0,
        "module_size_avg" : 0,
        "module_size_std" :0
    }
    if len(modules_summary) > 0:
        modules_summary = pd.DataFrame(modules_summary)
        disease_genes_extracted = float(len(set(module_genes).intersection(expected_genes)))
        disease_data["TP"] = disease_genes_extracted
        disease_data["TP/(TP+TN)_(_precision_)"] = disease_genes_extracted / len(module_genes)
        disease_data["TP/(TP+FN)_(_recall_)"] = disease_genes_extracted / len(expected_genes)
        if (disease_data["TP/(TP+TN)_(_precision_)"] + disease_data["TP/(TP+FN)_(_recall_)"]) == 0:
            disease_data["F1"] = 0
        else:
            disease_data["F1"] = 2 * ((disease_data["TP/(TP+TN)_(_precision_)"] * disease_data["TP/(TP+FN)_(_recall_)"]) /
                                  (disease_data["TP/(TP+TN)_(_precision_)"] + disease_data["TP/(TP+FN)_(_recall_)"]))

        disease_data["module_size_avg"] = modules_summary[SH_NUM_GENES].mean()
        disease_data["module_size_std"] = modules_summary[SH_NUM_GENES].std()


    pd.DataFrame([disease_data]).to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, dataset_name, algo_name, "{}_disease.tsv".format(report_file_name)),sep="\t", index=False)


def general_algo_report(algo_name, all_hg_reports, module_genes, modules_summary, report_file_name, total_hg_report, dataset_name):
    data = {}
    if len(modules_summary) > 0 :
        df_summary = pd.DataFrame(modules_summary)
        data = {"num_of_modules": df_summary.index.size,
                "module_size_avg": df_summary[SH_NUM_GENES].mean(),
                "module_size_std": df_summary[SH_NUM_GENES].std(),
                "total_num_genes": len(module_genes)
                }

    df = pd.DataFrame()
    if len(data) >0:
        df = pd.DataFrame([data])

    df.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, dataset_name, algo_name,
                     "{}_general.tsv".format(report_file_name)), sep="\t", index=False)



# def output_modules(output_file_name, modules, score_file_name, output_base_dir=""):
#     output_data = create_modules_output(modules, score_file_name)
#     file(output_file_name, 'w+').write(output_base_dir + "\n")
#     json.dump(output_data, file(output_file_name, 'a+'))
#     sys.stdout.write(output_file_name)

def reduce_to_dict(x,y):
    if y["id"] in x:
        x[y["id"]]["modules"] = x[y["id"]]["modules"] + y["modules"]
    else:
        x[y["id"]]=y
    return x

def merge_two_dicts(x, y):

    z = x.copy()
    z.update(y)
    return z

def create_modules_output(G_modules, score_file_name):
    scores=None
    if score_file_name is not None:
       print("score_file_name: {}".format(score_file_name))
       print(pd.read_csv(score_file_name,sep="\t").columns)
       scores = pd.read_csv(score_file_name,sep="\t").set_index("id")

       if constants.IS_PVAL_SCORES:
            scores["score"] = scores["pval"].apply(lambda x: -np.log10(x))

    zero_scores = [ {"score" : 0, "id" : gene} for G_module in G_modules for gene in G_module.nodes if scores is None or gene not in scores.index]
    if len(zero_scores) !=0:
        zero_scores = pd.DataFrame(zero_scores).set_index("id")
        zero_scores=zero_scores[~zero_scores.index.duplicated(keep='first')]
        scores = pd.concat([scores, zero_scores],axis=0)
    return [merge_two_dicts({"id" : k}, v) for k,v in reduce(reduce_to_dict, [{"eid": gene, "modules": [i], "id": gene, "gene_symbol": e2g_convertor([gene])[0], "score" : float(scores.loc[gene,"score"])} for i, G_module in enumerate(G_modules) for gene in G_module],\
            {}).items()]

def draw_network(G_modules, score_file_name, network_file_name):
    output = [{"data" : x, "label" : x["eid"], "selected" : True } for x in create_modules_output(G_modules, score_file_name)]
    # active_genes = [y for x in G_modules for y in x.nodes]
    active_edges = [y for x in G_modules for y in x.edges]
    # active_edges = [[x.iloc[0], x.iloc[2]] for i, x in pd.read_csv(network_file_name, sep="\t").iterrows() if x.iloc[0] in active_genes and x.iloc[2] in active_genes]
    additional_edges = [] # [[x.iloc[0], x.iloc[2]] for i, x in pd.read_csv(network_file_name, sep="\t").iterrows() if not (x.iloc[0] in active_genes and x.iloc[2] in active_genes) and (x.iloc[0] in active_genes or x.iloc[2] in active_genes)]
    # additional_nodes = [] # [y for x in (active_edges + additional_edges) for y in x if y if y not in active_genes]
    additional_nodes = [] # list(set(additional_nodes))

    return output + [{"data" : {"id" : x, "eid" : x, "modules" : []}, "label" : ""} for x in additional_nodes] + [{"data": {"id" : x[0]+"_"+x[1], "source":x[0], "target":x[1]}, "label" : ""} for x in additional_edges] + [{"data": {"id" : x[0]+"_"+x[1], "source":x[0], "target":x[1]}, "label" : "-"} for x in active_edges]



def generate_report_from_template(cy, output_base_dir, output_file_name):

    len([x for x in cy if not "source" in x["data"] and len(x["data"]["modules"])>0])
    report_file_name=format_script(os.path.join(os.path.dirname(os.path.abspath(__file__)),'../data', "graph.html"), NUM_OF_GENES=len([x for x in cy if not "source" in x["data"] and len(x["data"]["modules"])>0]), HG_REPORT=[], MODULES_SUMMARY=[], DISEASE_GENES=[], DATA=json.dumps(cy))

    shutil.move(report_file_name,
                os.path.join(output_base_dir, "module_{}.html".format(output_file_name)))
    return "module_{}.html".format(output_file_name)


def visualize_modules(dataset_name, G_modules, score_file_name, network_file_name, output_base_dir):
    print("visualizing modules...")
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)

    manager=multiprocessing.Manager()
    modules_summary = manager.list()

    params=[]
    for i, G_module in enumerate(G_modules):
        params.append([i, G_module, score_file_name, network_file_name, dataset_name, modules_summary, output_base_dir])
    p=multiprocessing.Pool(constants.N_OF_THREADS)
    p.map(module_report, params)
    p.close()
    # [module_report(p) for p in params]

def module_report(params):
    module_index, G_module, score_file_name, network_file_name, dataset_name, modules_summary, output_base_dir=params
    print("visualize module {} for dataset {}".format(module_index, dataset_name))

    modules_summary_row = {SH_MODULE_NAME: module_index, SH_NUM_GENES: len(G_module.nodes)}
    cy = draw_network([G_module], score_file_name, network_file_name)

    generate_report_from_template(cy, output_base_dir, str(module_index))
    if modules_summary is not None:
        modules_summary.append(modules_summary_row)
    return modules_summary


