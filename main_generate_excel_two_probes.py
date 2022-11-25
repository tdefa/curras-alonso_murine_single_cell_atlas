




import argparse
import os

from compute_spatial_state import generate_exels_two_probes

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='test')

    parser.add_argument("--path_save",
                        type=str,
                        default="/media/tom/250822/czi_folder/exels/",
                        help='path to save the exels')

    parser.add_argument("--path_folders_czi",
                        type=str,
                        default="/media/tom/250822/czi_folder/",
                        help='Path to the  folder with all the experiment folder contiaing the czi file')

    parser.add_argument("--list_folder", nargs="+", default=["00_Macrophages/"],
                        help='List of folders to analyse in the czi main folder')


    parser.add_argument('--list_probes1',
                        type=str, nargs='+', action='append', default=[['Chil3']],
                        help="")

    parser.add_argument('--list_probes2', type=str, nargs='+', action='append', default=[
        ['Mki67']])

    parser.add_argument("--dico_stats_name",
                        type=str,
                        default="finaldico_241122.npy",
                        help='name of the exels file')


    parser.add_argument("--compute_nuclei_size",
                        type=bool,
                        default=True,
                        help='if true it computes also the nuclei size (can take some time)')

    parser.add_argument("--scale_z",
                        type=float,
                        default=300,
                        help='if true it computes also the nuclei size')

    parser.add_argument("--scale_xy",
                        type=float,
                        default=103,
                        help='if true it computes also the nuclei size')

    parser.add_argument("--port", default=39948)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')

    args = parser.parse_args()
    print(args)

    assert len(args.list_probes1) == len(args.list_probes2)

    try:
        os.mkdir(args.path_save)
    except Exception as e:
        print(e)
    for probes_index in range(len(args.list_probes1)):
        probes1 = args.list_probes1[probes_index]
        probes2 = args.list_probes2[probes_index]
        generate_exels_two_probes(list_folder=args.list_folder,
                                       gene_1=probes1,
                                       gene_2=probes2,
                                       path_to_take=args.path_folders_czi,
                                       path_save=args.path_save,
                                       dico_stat_name=args.dico_stats_name,
                                       compute_nuclei_size=args.compute_nuclei_size,
                                       scale_z=args.scale_z,
                                       scale_xy=args.scale_xy,
                                       )


