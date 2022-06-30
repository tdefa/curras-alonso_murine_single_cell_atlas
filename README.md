

# All the usefull function and the pipeline is in main_cluster_split

# to run it:

##### activate your conda env ex: 

#### execute the following command in the console in the code directory it will compute all the mapping and images
python main_cluster_split.py --path_to_czi_folder /media/tom/Elements1/to_take/test_pipeline/ --list_folder test1/ test2/ --new_probe Pdgfratest 35 0.42 --new_probe Pdgfratest2 40 0.62  --dico_name_save test_for_hugo --prepare_czi 1 --segmentation 1 --spot_detection 1 --classify 1 --save_plot 1

### generate the exels files for one cell type use the following command

python compute_spatial_state --path_save /media/tom/Elements1/to_take/test_pipeline/exels/ --path_folder_of_folders /media/tom/Elements1/to_take/test_pipeline/ --list_folder test1 test2 --list_probes Lamp3 --list_probes 'Cap' 'aCap' 'CEC' 'acap' --dico_stat_name finaltest_for_hugo.npy


# Argument of the commend line 

optional arguments:
  -h, --help            show this help message and exit
  -ptz PATH_TO_CZI_FOLDER, --path_to_czi_folder PATH_TO_CZI_FOLDER
                        path to the folder containing the czi
                        
  --list_folder LIST_FOLDER [LIST_FOLDER ...]
                        list of folders in the czi folders to analyse
                        
  --new_probe NEW_PROBE [NEW_PROBE ...]
                        command to add new probes or change parameters of existing one to add it do --new_probe p1 epsi overlapping
                        --new_probe p2 20 0.3 where 'epsi' is the parameter of the dbscan 'overlapping' is the percent of overlap to make
                        a cell positive to a probe
                        
  --manual_threshold_cy3 MANUAL_THRESHOLD_CY3
                        write a json like the : {"02_NI1230_Lamp3-Cy5_Pdgfra-Cy3_08.tiff": 8, "01_IR5M1236_Lamp3-Cy5_Pdgfra-Cy5_04.tiff":
                        7} to set manually the rna spot detection threshold
                        
                        
  --manual_threshold_cy5 MANUAL_THRESHOLD_CY5
  
  
  
  -dns DICO_NAME_SAVE, --dico_name_save DICO_NAME_SAVE
                        additional name in the save result

  
  
   --prepare_czi    do : prepare_czi to tiff 1 , do not 0
  -sg SEGMENTATION, --segmentation SEGMENTATION
                        do segmentation 1 , do not 0
  --spot_detection SPOT_DETECTION
                        do spots detection 1 , do not 0
  --classify CLASSIFY   do classification / cell type mapping 1 , do not 0
  --save_plot SAVE_PLOT do save plot 1 , do not 0




###################

helping function 
