import src.constants as constants
import os
g2e_dict = None
e2g_dict = None
dict_type=""

def load_gene_dictionary(gene_list_file_name, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path is None:
        gene_list_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../data",gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines


def get_g2e_dictionary(cur_dict_type='ENSG'):

    if cur_dict_type=='ENSMUSG':
        lines_dict = load_gene_dictionary(constants.ENSMUSG_TO_GENE_SYMBOLS)
    elif cur_dict_type=='ENSG':
        lines_dict = load_gene_dictionary(constants.ENSG_TO_GENE_SYMBOLS)
    else:
        raise "unknown gene identifiers: {}".format(cur_dict_type)

    gene_symbols2ensembl = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
        gene_symbols2ensembl[splited_line[1]] = splited_line[0][:limit]
    return gene_symbols2ensembl

def get_e2g_dictionary(cur_dict_type):
    if cur_dict_type=='ENSMUSG':
        lines_dict = load_gene_dictionary(constants.ENSMUSG_TO_GENE_SYMBOLS)
    elif cur_dict_type=='ENSG':
        lines_dict = load_gene_dictionary(constants.ENSG_TO_GENE_SYMBOLS)
    else:
        raise "unknown gene identifiers: {}".format(cur_dict_type)

    ensembl2gene_symbols = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
        ensembl2gene_symbols[splited_line[0][:limit]] = splited_line[1]

    global dict_type
    dict_type=cur_dict_type
    return ensembl2gene_symbols


def e2g_convertor(e_ids):
    if type(e_ids) is str:
        e_ids=[e_ids]

    global g2e_dict
    global dict_type
    cur_dict_type=e_ids[0][:e_ids[0].index('0')]
    if g2e_dict is None or dict_type!=cur_dict_type:
        e2g_dict = get_e2g_dictionary("ENSG") # cur_dict_type)

    results = []
    for cur in e_ids:
        if cur.split(".")[0] in e2g_dict:
            results.append(e2g_dict[cur.split(".")[0]])
        else:
            results.append(cur.split(".")[0])
    return results
