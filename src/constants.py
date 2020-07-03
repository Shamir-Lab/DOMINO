import json
import os

from matplotlib.lines import Line2D
import numpy as np

import matplotlib as mpl
font = {'size'   : 22}
mpl.rc('font', **font)
mpl.rc('xtick', labelsize=30)    # fontsize of the tick labels
mpl.rc('ytick', labelsize=30)
mpl.rc('axes', labelsize=22)
mpl.rc('legend', fontsize=20)




dir_path = os.path.dirname(os.path.realpath(__file__))
PATH_TO_CONF = "env/config/conf.json"

REPO_DIR = os.path.dirname(os.path.realpath(__file__))
SH_DIR = os.path.join(REPO_DIR, "sh","scripts")


LABEL_ID = "sample_type.samples"
PRIMARY_TUMOR = "Primary Tumor"
METASTATIC = "Metastatic"

LABELS_NORMAL = "labels_normal"
LABELS_SHUFFLE = "labels_shuffle"
LABELS_RANDOM = "labels_random"
LABELS_ALTERNATED = "labels_alternated"
LABELS_INVERTED = "labels_inverted"

ENSG_TO_GENE_SYMBOLS = "ensg2gene_symbol.txt"
ENSMUSG_TO_GENE_SYMBOLS = "ensmusg2gene_symbol.txt"
ENSEMBL_TO_ENTREZ = "ensembl2entrez.txt"

GO_OBO_URL = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
GO_ASSOCIATION_GENE2GEO_URL = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'
GO_FILE_NAME = 'go_bp.obo' #'go-basic.obo'
GO_ASSOCIATION_FILE_NAME = "gene2go"

