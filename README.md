# DOMINO

DOMINO: Discovery of Modules In Networks using Omic.

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
    git clone https://github.com/hag007/domino.git
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
slicer --datasets </path/to/dataset1,/path/to/dataset2...> --output_dir <prefix of dataset>
```

To run DOMINO:
```
domino --active_genes_files </path/to/dataset1,/path/to/dataset2...> --network_file </path/to/network.sif> --slices_file <slices_file.txt> --output_folder </path/to/output_folder> [--sth <slices_threshold> --mth <putative_modules_threshold>]
```

The common command line options are:

`-a/--active_genes_files`: list of files of active genes. gene ids are Enseble id, separated by new line.

`-n/--network_file`: path to network file (sif format).

`-s/--slices_file`: path to slices file (i.e. the output of "slicer" script).


## Main output files

`output_folder/active_gene_file_name.out`: list of final modules

## Advanced usage

`-sth/--slices_threshold`: threshold for considering a slice as relevant

`-mth/--module_threshold`: threshold for considering a putative module as final module.
