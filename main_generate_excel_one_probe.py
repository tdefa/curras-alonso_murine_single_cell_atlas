

import argparse
from compute_spatial_state import generate_exels_one_cell

if __name__ == '__main__':
    #%% one cell
    list_folder = [
        "test1/"
    ]

    parser = argparse.ArgumentParser(description='test')


    parser.add_argument("--path_save",
                        type=str,
                        default="/media/tom/Elements1/to_take/test_pipeline/exels/",
                        help='path to save the exels')

    parser.add_argument("--path_folder_of_folders",
                        type=str,
                        default="/media/tom/Elements1/to_take/test_pipeline/",
                        help='path to the  folder with all the experiment')


    parser.add_argument("--list_folder", nargs="+", default=list_folder,  help=' list of folders in the czi folders to analyse') #


    parser.add_argument('--list_probes', type=str, nargs='+', action='append', default=[ ['Lamp3'],  ['Pecam1'],  ['Ptprb'],['Hhip'], ["Rtkn2"],
                 ['Apln'], ['Chil3'],  ['Fibin'], ['C3ar1'],
                  ['Mki67'],  ['Cap', 'aCap', 'CEC', 'acap'],])

    parser.add_argument("--dico_stat_name",
                        type=str,
                        default="finaltest_for_hugo.npy",
                        help='path to save the exels')


    parser.add_argument("--port", default=39949)
    parser.add_argument("--mode", default='client')

    args = parser.parse_args()
    print(args)

    for probes in args.list_probes:
        generate_exels_one_cell(list_folder =args.list_folder,
                            gene_smfish =probes,
                            path_to_take = args.path_folder_of_folders,
                            path_save = args.path_save,
                            dico_stat_name = args.dico_stat_name,
                            compute_nuclei_size = False,
                            nuclei_shape_index = False,
                            compute_neighbor = False)

