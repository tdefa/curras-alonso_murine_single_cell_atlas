

import argparse

from compute_spatial_state import generate_exels_one_cell

if __name__ == '__main__':
    #%% one cell

    parser = argparse.ArgumentParser(description='test')
    parser.add_argument("--path_save",
                        type=str,
                        default="/media/tom/250822/czi_folder/exels/",
                        help='Path to the folder where the exels containing cell type calling results will be saved')

    parser.add_argument("--path_folders_czi",
                        type=str,
                        default="/media/tom/250822/czi_folder/",
                        help='Path to the parent folder with all the experiment folder containg the czi file')
    parser.add_argument("--list_folder", nargs="+", default=["00_Macrophages/"],
                        help='List of folders to analyse in the czi main folder')
    #parser.add_argument("--list_probes", nargs="+", default=[],
    #                    help="Probes to analyse")
    parser.add_argument('--list_probes', type=str, nargs='+', action='append',
                        default=[],
                        help="Probes to analyse")

    parser.add_argument("--dict_name_save",
                        type=str,
                        default="finaldico_241122.npy",
                        help='name of the .npy file containing the analysis result')

    parser.add_argument("--compute_nuclei_size",
                        type=bool,
                        default=True,
                        help='if true it computes also the nuclei size (can take some time)')

    parser.add_argument("--scale_z",
                        type=float,
                        default=300,
                        help='pixel z size in nm')

    parser.add_argument("--scale_xy",
                        type=float,
                        default=103,
                        help='pixel xy size in nm')

    parser.add_argument("--port", default=39948)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')




    args = parser.parse_args()
    print(args)

    if args.path_save[-1] != '/':
        args.path_save[-1] += '/'
    if args.path_folders_czi[-1] != '/':
        args.path_folders_czi[-1] += '/'
    print(f'list probe {args.list_probes}')

    for probes in args.list_probes:
        generate_exels_one_cell(list_folder=args.list_folder,
                                gene_smfish=probes,
                                path_to_take=args.path_folders_czi,
                                path_save=args.path_save,
                                dico_stat_name=args.dict_name_save,
                                compute_nuclei_size=args.compute_nuclei_size,
                                scale_z=args.scale_z,
                                scale_xy=args.scale_xy
                                )


