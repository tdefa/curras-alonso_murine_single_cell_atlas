

python main_cluster.py --path_to_czi_folder ./demo_data/ \
--list_folder experiment1/  \
--dict_name_save analysis2022  \

python  main_generate_excel_one_probe.py \
--path_save ./demo_data/excel/ \
--path_folders_czi ./demo_data/ \
--list_folder experiment1/ \
--list_probes 'Chil3'   \
--list_probes 'Mki67' \
--dict_name_save analysis2022.npy \



python  main_generate_excel_two_probes.py \
--path_save ./demo_data/excel/ \
--path_folders_czi ./demo_data/ \
--list_folder experiment1/ \
--list_probes1 'Chil3' \
--list_probes2 'Mki67' \
--dict_name_save analysis2022.npy \
