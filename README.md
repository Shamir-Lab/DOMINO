# DOMINO

DOMINO: Discovery of Modules In Networks using Omic.

DOMINO is a network-based module discovery (NBMD) algorithm.  It recieves a gene network and nodes’ activity scores as input and report sub-networks (modules) that are putatively biologically meaningful in the context of the activity data.


In extensive evaluation conducted on gene expression and genome-wide association study data we discovered that NBMD algorithms tended to over-reporting of enrichment: GO terms enriched in the modules on real data were often also enriched when the algorithms were run on randomly permuted activity scores.

In constrast, modules retrieved by DOMINO had high rate of empirically validated GO terms.

A preprint version of the study is available at https://www.biorxiv.org/content/10.1101/2020.03.10.984963v3.

- [Installation](#installation)
- [Basic Usage](#basic-usage)
- [Main output files](#main-output-files)
- [Advanced usage](#advanced-usage)

## Installation

### From sources
Download the sources and install according to the following:

DOMINO is written in Python3. The necessary libraries will all be installed by the `setup.py` script.
We recommend using a virtual environment. For example, in Linux, before running `setup.py`:
```
python -m venv domino-env
source domino-env/bin/activate
```
To install, download and run setup.py:
```
    git clone https://github.com/Shamir-Lab/DOMINO.git
    cd DOMINO
    python setup.py install
```
It is possible to install as a user without root permissions:
```
python setup.py install --user
```

## Basic Usage

To run preprocessing step 0 (partitioning network according Newman-Girvan algorithm):
```
slicer --network_file </path/to/network.sif> --output_dir <prefix of dataset>
```

To run DOMINO:
```
domino --active_genes_files </path/to/dataset1,/path/to/dataset2...> --network_file </path/to/network.sif> --slices_file <slices_file.txt> --output_folder </path/to/output_folder> [-sth <slices_threshold> -mth <putative_modules_threshold>]
```

The common command line options are:

`-a/--active_genes_files`: list of files of active genes. gene ids are Enseble id, separated by new line.

`-n/--network_file`: path to network file (sif format).

`-s/--slices_file`: path to slices file (i.e. the output of "slicer" script).


## Main output files

`output_folder/active_gene_file_name/modules.out`: list of final modules
`output_folder/active_gene_file_name/module_i.html`: visualization of the i'th module

## Advanced usage

`-sth/--slices_threshold`: threshold for considering a slice as relevant

`-mth/--module_threshold`: threshold for considering a putative module as final module.
