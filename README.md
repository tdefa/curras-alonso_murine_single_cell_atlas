



#### Module to infer the cell type using smFISH data with two channel of two different marker genes. All the usefull functions and the pipeline is in main_cluster_split


#### 1) create your conda from the file ENV.txt: 
conda list --explicit>ENV.txt

#### 2) execute the following command in the console in the code directory it will compute all the mapping and images. Replace by the parameter by your configuration.  run  python main_cluster_split.py --help for more information

python main_cluster_split.py --path_to_czi_folder /media/tom/Elements1/to_take/test_pipeline/ --list_folder test1/ test2/ --new_probe Pdgfratest 35 0.42 --new_probe Pdgfratest2 40 0.62  --dico_name_save test_for_hugo    --manual_threshold_cy3 {"02_NI1230_Lamp3-Cy5_Pdgfra-Cy3_08.tiff": 8, "01_IR5M1236_Lamp3-Cy5_Pdgfra-Cy5_04.tiff": 7} --prepare_czi 1 --segmentation 1 --spot_detection 1 --classify 1 --save_plot 1

### 3) generate the exels files for one cell type use the following command

python compute_spatial_state --path_save /media/tom/Elements1/to_take/test_pipeline/exels/ --path_folder_of_folders /media/tom/Elements1/to_take/test_pipeline/ --list_folder test1 test2 --list_probes Lamp3 --list_probes 'Cap' 'aCap' 'CEC' 'acap' --dico_stat_name finaltest_for_hugo.npy


