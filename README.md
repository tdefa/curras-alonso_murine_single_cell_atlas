



#### Module to infer the cell type using smFISH data with two channels of two different marker genes. All the usefull functions and the pipeline is in main_cluster_split


#### 1) create your conda from the file ENV.txt: 
conda list --explicit>ENV.txt

#### 2) execute the following command in the console in the code directory it will compute all the mapping and images. Replace by the parameter by your configuration.  run  python main_cluster_split.py --help for more information

python main_cluster_split.py --path_to_czi_folder [path to the folder containing folders of czi of the defferent experience] \
--list_folder experience1/ experience1/ \ #[name of folder in the czi_folder containing the czi image file you want to analyse]
--new_probe Pdgfratest 35 0.42 #[probe name, espi parameter in dbscan, minimal overlapping to make a nucleus positive]
--new_probe Pdgfratest2 40 0.62 
--dico_name_save analyisi2022 #[key word use to name files storing result of the analysis]  


### 3) generate the exels files for one cell type use the following command

python compute_spatial_state --path_save /media/tom/Elements1/to_take/test_pipeline/exels/ \
--path_folder_of_folders /media/tom/Elements1/to_take/test_pipeline/  \
--list_folder test1 test2  \ 
--list_probes Lamp3 \ 
--list_probes 'Cap' 'aCap' 'CEC' 'acap' \ 
--dico_stat_name finaltest_for_hugo.npy


