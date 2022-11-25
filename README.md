



#### Module to infer the cell type using smFISH data with two channels of two different marker genes. All the usefull functions and the pipeline is in main_cluster_split


#### 1) create your conda from the file ENV.txt: 
conda list --explicit>ENV.txt

#### 2) execute the following command in the console in the code directory it will compute all the mapping and images. Replace by the parameter by your configuration.  run  python main_cluster_split.py --help for more information

python main_cluster_split.py --path_to_czi_folder /media/tom/250822/czi_folder/ #[path to the folder containing folders of czi of the defferent experience] \
--list_folder experience1/ experience1/ \ #[name of folder in the czi_folder containing the czi image file you want to analyse] <br />
--new_probe Pdgfratest 35 0.42 \ #[probe name, espi parameter in dbscan, minimal overlapping to make a nucleus positive] <br />
--new_probe Pdgfratest2 40 0.62 \ <br />
--dico_name_save analyisi2022 #[key word use to name files storing result of the analysis]  <br />


### 3) generate the exels files for one cell type use the following command

python  main_generate_excel_one_probe.py \ <br />
--path_save /media/tom/Elements1/to_take/test_pipeline/exels/ \ #[Path to the folder where the exels containing cell type calling result will be saved] <br />
--path_folders_czi /media/tom/250822/czi_folder/  #[Path to the  folder with all the experiment folder contiaing the czi file] \ <br />
--list_folder test1/ test2/ #[List of folders to analyse in the czi main folder]
--list_probes Lamp3 \ #[Probes to analyse, add space to separate probes that have many names] <br />
--list_probes 'Cap' 'aCap' 'CEC' 'acap' \ #[Probes to analyse, add space to separate probes that have many names (like here)] <br />
--dico_stat_name .npy \ #[Name of the dictionary storing cell type calling result from main_cluster.py] <br />


### 4) generate the exels files to compare size and proportion of two diffenrent cell type.



python  main_generate_excel_one_probe.py \ <br />
--path_save /media/tom/Elements1/to_take/test_pipeline/exels/ \ #[Path to the folder where the exels containing cell type calling result will be saved] <br />
--path_folders_czi /media/tom/250822/czi_folder/  #[Path to the  folder with all the experiment folder contiaing the czi file] \ <br />
--list_folder test1/ test2/ #[List of folders to analyse in the czi main folder]
--list_probes1 Lamp3 \ #[First set of probes to analyse] <br />
--list_probes2 'Cap' 'aCap' 'CEC' 'acap' \ #[second set of probes that will be compared the the first set, add space to separate probes that have many names (like here)]  <br />
--list_probes2 Mki67
--dico_stat_name finaldico_241122.npy \ #[Name of the dictionary storing cell type calling result from main_cluster.py] <br />


### use python [name of the command] --help for a more extensive documentation



