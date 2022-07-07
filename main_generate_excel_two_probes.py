




import os
import argparse
from compute_spatial_state import generate_exels_cell_state_type



if __name__ == '__main__':


    list_folder = [
        #"220218_IR3M17gy/",
        #"220317_IR3M10gy_2346_A/",
        #"220318_IR3M10gy_2346_C/",
        #"220325_IR3M17gy_2351_B/",
        #"220425_IR3M10gy_2345_B/",
        #"220502_fibro_gcap_B/",
        #"220304_IR3M17gy_prolicence/",
        #"220317_IR3M10gy_2346_B/",
        #"220324_IR3M17gy_2351_A/",
        #"220422_IR3M10gy_2345_A/",
        #"220429_fibro_gcap_A/",
        #"220503_fibro_gcap_C/",
        #"220505_fibro_gcap_D/",
        #"220506_fibro_gcap_E/",
        "test1/"
    ]
    parser = argparse.ArgumentParser(description='test')

    parser.add_argument("--path_save",
                        type=str,
                        default="/media/tom/Elements1/to_take/test_pipeline/exels_two_probes/",
                        help='path to save the exels')

    parser.add_argument("--path_folder_of_folders",
                        type=str,
                        default="/media/tom/Elements1/to_take/test_pipeline/",
                        help='path to the  folder with all the experiment')

    parser.add_argument('--list_probes1', type=str, nargs='+', action='append', default=[
        ['Cap', 'aCap', 'CEC',  'acap'],['Cap', 'aCap', 'CEC',  'acap'], ["Fibin"]],
                        help = "")

    parser.add_argument('--list_probes2', type=str, nargs='+', action='append', default=[
        ['Mki67'], ['Cap', 'aCap', 'CEC', 'acap'], ['Serpine1']])

    parser.add_argument("--dico_stats_name",
                        type=str,
                        default="finaltest_for_hugo.npy",
                        help='path to save the exels')



    parser.add_argument("--port", default=39949)
    parser.add_argument("--mode", default='client')
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

        generate_exels_cell_state_type(list_folder = list_folder,
                                       gene_1 = probes1,
                                       gene_2 = probes2 ,
                                       path_to_take = args.path_folder_of_folders,
                                       path_save = args.path_save,
                                       dico_stat_name = args.dico_stats_name,
                                       scale_z=300,
                                       scale_xy=103,
                                       compute_morpho_features=False)


