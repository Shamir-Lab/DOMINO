import argparse
import os
from src.core.domino import main as domino_main
from src.core.preprocess_slices import create_slices
from src.utils.visualize_modules import visualize_modules

def main_domino():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('-a', '--active_genes_files', dest='active_genes_files', help='/path/to/active_genes_files_1,/path/to/active_genes_files_2', default="/media/hag007/Data/DOMINO_amit/brca.txt")
    parser.add_argument('-o', '--output_folder', dest='output_folder', help='/path/to/output', default="/media/hag007/Data/domino_test/out_test/")
    parser.add_argument('-n', '--network_file', dest='network_file', help='/path/to/network file', default="/media/hag007/Data/DOMINO_amit/network/dip.sif")
    parser.add_argument('-s', '--slices_file', dest='slices_file', help='/path/to/slices file', default="/media/hag007/Data/DOMINO_amit/network/dip_sliced.txt")
    parser.add_argument('-sth', '--slice_threshold', dest='slice_threshold', default="0.3", help='threshold of slices')
    parser.add_argument('-mth', '--module_threshold', dest='module_threshold', default="0.05", help='threshold of putative modules')


    args = parser.parse_args()
    active_genes_files = args.active_genes_files.split(",")
    output_folder = args.output_folder
    network_file = args.network_file
    slices_file = args.slices_file
    slice_threshold = float(args.slice_threshold)
    module_threshold = float(args.module_threshold)

    for cur_ag in active_genes_files:
        final_modules=domino_main(active_genes_file=cur_ag, network_file=network_file, slices_file=slices_file, slice_threshold=slice_threshold, module_threshold=module_threshold)
        activity_name=os.path.splitext(os.path.split(cur_ag)[-1])[0]
        report_folder=os.path.join(output_folder,activity_name)
        try:
            os.makedirs(report_folder)
        except:
            pass

        out_file=os.path.join(report_folder, "modules.out")
        open(out_file, 'w+').write("\n".join(['[%s]' % ', '.join(m) for m in final_modules]))
        print(f'{len(final_modules)} final modules are reported at {out_file}')

        visualize_modules(os.path.splitext(cur_ag.split('/')[-1])[0], final_modules, None, network_file, report_folder)

def main_slicer():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('-n', '--network_file', dest='network_file', help='/path/to/network', default="/media/hag007/Data/bnet/networks/dip_reduced.sif")
    parser.add_argument('-o', '--output_file', dest='output_file', default="/media/hag007/Data/bnet/networks/dip_sliced.txt", help='/path/to/output')


    args = parser.parse_args()
    network_file = args.network_file
    output_file = args.output_file
    create_slices(network_file, output_file)




if __name__=="__main__":
    # main_slicer()
    main_domino()