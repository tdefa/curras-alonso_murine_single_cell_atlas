#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from utils_ext.cellpose_utilis import stitch3D
import os
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
import tifffile
import numpy as np
import cellpose
from cellpose import models, io, plot

import argparse
def segment_nuclei(path_to_dapi, path_to_mask_dapi, dico_param, model, save = True):
    onlyfiles = [f for f in listdir(path_to_dapi) if isfile(join(path_to_dapi, f))]
    print(onlyfiles)
    for f in onlyfiles:
        print(f)
        img = tifffile.imread(path_to_dapi + f)
        print(img.shape)
        if dico_param["gpu"] and (img.shape[-1] > 6000 or img.shape[-2] > 6000):
            continue

        if dico_param["mip"] is True and len(img.shape) == 3:
            img = np.amax(img, 0)
       # elif dico_param["mip"]
        else:
            if len(img.shape) == 3:
                img = img.reshape(img.shape[0], 1, img.shape[1], img.shape[2])
                print(f'image dapi shape after reshape {img.shape}')
                img = list(img)

        masks, flows, styles, diams = model.eval(img, diameter=dico_param["diameter"],
                                                 channels=[0,0],
                                                flow_threshold=dico_param["flow_threshold"], do_3D= dico_param["do_3D"],
                                                stitch_threshold = 0)
        try:
            # I stich with my own custum funvtion to avoid  "zero-size array to reduction operation maximum which has no identity"
            masks = stitch3D(masks, dico_param["stitch_threshold"])

        except Exception as e:
                    print(e)
                    print("the file %s lead to an error" % f)
                    print()
                    break
        masks = np.array(masks)
        if len(masks.shape) < 3:
            plt.imshow(masks * 3)
            plt.show()
            plt.imshow(img)
            plt.show()
        if save:
            tifffile.imwrite(path_to_mask_dapi + "dapi_mask" + f[:-4] + "tiff", data=masks, dtype=masks.dtype)
            np.save(path_to_mask_dapi + "dico_param.npy", dico_param)

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='test')
    ### path to dapi tiff datset
    parser.add_argument('-pi',"--path_input", type=str, default= "/media/tom/Transcend/microscope-comparison/widefield/dapi/" , help='')
    parser.add_argument('-po',"--path_output", type=str, default= "/media/tom/Transcend/microscope-comparison/widefield/dapi_mask/" , help='')

   # parser.add_argument()

    ###cellpose arg
    parser.add_argument('-d',"--diameter", type=float, default= None, help='')
    parser.add_argument('-ft',"--flow_threshold", type=float, default=0.8, help='')
    parser.add_argument('-d3',"--do_3D", type=bool, default= False, help='')
    parser.add_argument('-m',"--mip", type=bool, default=False, help='')
    parser.add_argument('-st',"--stitch_threshold", type=float, default= 0.4, help='')
    #parser.add_argument('-mode',"--mode", type=str, default='nuclei', help='')
    parser.add_argument("--port", default=39949)
    parser.add_argument("--mode", default='client')

    args = parser.parse_args()
    
    folder_name =  "predicted_mask" +"/"

    #if not os.path.exists(args.path_output + folder_name):
     #    os.mkdir(args.path_output + folder_name)
    model = models.Cellpose(gpu=False, torch =True, model_type='nuclei')
    # ##parameter
    dico_param = {}
    dico_param["diameter"]  = args.diameter
    dico_param["flow_threshold"] = args.flow_threshold
    dico_param["do_3D"] = args.do_3D
    dico_param["mip"] = args.mip
    dico_param["projected_focused"] = False
    dico_param["stitch_threshold"] = args.stitch_threshold

    r = segment_nuclei(args.path_input ,
                   args.path_output ,
                   dico_param, model)



#%%
    list_folder = [
     "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Capillary_EC/", #
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Large_Vessels/", #reaploaded
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Macrophages/", #reaploaded
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200908_CEC/", #reaploaded
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200908_fibrosis/", #decalage
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/201030_fridyay/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/201127_AM_fibro/", ##reaploadedss
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210205_Prolicence/aCap_prolif/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210205_Prolicence/aCap_senes/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210219_myo_fibros_y_macrophages/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210412_repeat_fibro/IR5M/", #decalage
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210412_repeat_fibro/NI/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210413_rep2/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210425_angiogenesis/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210426_repeat3/",
    
    ]
    
    image_to_recompute =["01_IR5M1236_Pdgfra-Cy5_Hhip-Cy3_mid_02.tiff"]
        
    bug_image = []    
                       
                         
    for folder in list_folder:
        path_output_segmentaton = folder + "/tiff_data/predicted_mask_dapi/"
        onlyfiles = [f for f in listdir(path_output_segmentaton) if isfile(join(path_output_segmentaton, f)) and f[-1] == "f" ]
        onlyfiles = [onlyfiles[i][14:] for i in range(len(onlyfiles))]
        for f in onlyfiles:

            print(f)
            img = tifffile.imread(folder+"/tiff_data/dapi/dapi_"+f  )
            img = list(img.reshape(54,1, 1040, 1388))
            model = models.Cellpose(gpu=True, torch =True, model_type="nuclei")
            
            tifffile.imread = path_output_segmentaton + "dapi_maskdapi_" + f[:-4] + "tiff"
            i = -1
            while not np.sum(mask[i]>0) > 0:
                i -=1
            if not max(np.unique(masks[i])) + 1 == len(np.unique(masks[i])):
                continue
            bug_image.append(folder, f)
            masks, flows, styles, diams = model.eval(img, diameter=None, channels = [0,0],
                                                                do_3D= False,
                                                                stitch_threshold = 0.4)
                                                                #anisotropy = 3)
                                                                
            if  max(np.unique(masks[i])) + 1 == len(np.unique(masks[i])):
                masks = stitch3D(masks, 0.4)
            else: 
                print("error")
                break
            
            img = tifffile.imread(folder+"/tiff_data/dapi/dapi_"+f  )

            plt.imshow(np.amax(img,0))
            plt.show()
            
            masks = np.array(masks)
            plt.imshow(np.amax(masks, 0))
            plt.title(f)
            plt.show()
            tifffile.imwrite(path_output_segmentaton + "dapi_maskdapi_" + f[:-4] + "tiff", data=masks, dtype=masks.dtype)
        
            print("save")
            print(f)

#%%
    list_folder = [
     "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Capillary_EC/", #
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Large_Vessels/", #reaploaded
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Macrophages/", #reaploaded
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200908_CEC/", #reaploaded
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200908_fibrosis/", #decalage
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/201030_fridyay/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/201127_AM_fibro/", ##reaploadedss
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210205_Prolicence/aCap_prolif/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210205_Prolicence/aCap_senes/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210219_myo_fibros_y_macrophages/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210412_repeat_fibro/IR5M/", #decalage
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210412_repeat_fibro/NI/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210413_rep2/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210425_angiogenesis/",
    "/home/tom/Bureau/annotation/cell_type_annotation/to_take/210426_repeat3/",
    
    ]
#%%    
    list_folder = ["/home/tom/Bureau/annotation/cell_type_annotation/to_take/210428_IR5M1236_Lamp3-Cy5_Pdgfra-Cy3/"]
    
        
    bug_image = []    
                       
                         
    for folder in list_folder:
        path_output_segmentaton = folder + "/tiff_data/predicted_mask_dapi/"
        onlyfiles = [f for f in listdir(path_output_segmentaton) if isfile(join(path_output_segmentaton, f)) and f[-1] == "f" ]
        onlyfiles = [onlyfiles[i][14:] for i in range(len(onlyfiles))]
        for f in onlyfiles:

            
            masks = tifffile.imread(path_output_segmentaton + "dapi_maskdapi_" + f[:-4] + "tiff")
            i = -1
            while not np.sum(masks[i]>0) > 0:
                i -=1
            if not max(np.unique(masks[i])) + 1 == len(np.unique(masks[i])):
                continue
            print((folder, f))

            bug_image.append((folder, f))
            masks = stitch3D(masks, 0.4)
            tifffile.imwrite(path_output_segmentaton + "dapi_maskdapi_" + f[:-4] + "tiff", data=masks, dtype=masks.dtype)
#%%
        [ "10_IR5M_Ptprb-Cy3_Mki67-Cy5_05",
                         "10_IR5M_Ptprb-Cy3_Mki67-Cy5_06", 
                         "10_IR5M_Ptprb-Cy3_Mki67-Cy5_07",
                         "11_NI_Ptprb-Cy3_Serpine1-Cy5_01",
                         "11_NI_Ptprb-Cy3_Serpine1-Cy5_02", 
                         "11_NI_Ptprb-Cy3_Serpine1-Cy5_03",
                         "11_NI_Ptprb-Cy3_Serpine1-Cy5_04",
                         "11_NI_Ptprb-Cy3_Serpine1-Cy5_05",
                         "12_IR5M_Ptprb-Cy3_Serpine1-Cy5_01"
                         "12_IR5M_Ptprb-Cy3_Serpine1-Cy5_02",
                         "12_IR5M_Ptprb-Cy3_Serpine1-Cy5_05",
                         "02_IR5M_Chil3-Cy3_Mki67-Cy5_01",
                         "04_IR5M_Hhip-Cy3_Pdgfra-Cy5_002",
                         "IR1M_aCapCy3_Mki67Cy5_06", 
                         "Ctrl_aCapCy3_Mki67Cy5_07",
                         "IR4M_aCapCy3_Mki67Cy5_07",
                         "03_NI_Chil3-Cy3_Serpine1-Cy5_003",
                         "02_IR5M_Lamp3-Cy3_Pdgfra-Cy5_024",
                         "02_IR4M_Lamp3-Cy5_Pdgfra-Cy3_01",
                         "04_IR5M2201()_Pecam1-Cy5_Ptprb-Cy3_05",
                         "01_IR5M1236_Pdgfra-Cy5_Hhip-Cy3_mid_03",
                         '03_IR5M1249_Lamp3-Cy5_Pdgfra-Cy3_perif_07',
                         "04_IR5M1249_Pdgfra-Cy5_Hhip-Cy3_mid_05", 
                         "04_IR5M2201()_Pecam1-Cy5_Ptprb-Cy3_15",
                         "04_IR5M2201()_Pecam1-Cy5_Ptprb-Cy3_09",
                         "04_IR5M2201()_Pecam1-Cy5_Ptprb-Cy3_10",
                         "04_IR5M2201()_Pecam1-Cy5_Ptprb-Cy3_06",
                         "04_IR5M2201()_Pecam1-Cy5_Ptprb-Cy3_07",
                         "04_IR5M2201()_Pecam1-Cy5_Ptprb-Cy3_04",
                         "04_IR5M2201()_Pecam1-Cy5_Ptprb-Cy3_11",
                         "06_NI1230_Chil3-Cy5_C3ar1-Cy3_04",
                         "03_NI1225_Pdgfra-Cy5_Hhip-Cy3_04",
                         "03_NI1225_Pdgfra-Cy5_Hhip-Cy3_11",
                         "07_IR5M2330_Lamp3-Cy5_Pdgfra-Cy3_15",]

"""# Save
dictionary = {'hello':'world'}
np.save('my_file.npy', dictionary)

# Load
read_dictionary = np.load('/home/thomas/Bureau/phd/first_one/tiff_data/predicted_mask_dapi/dico_param.npy',allow_pickle='TRUE').item()
print(read_dictionary['hello']) # displays "world



fig = plt.figure(figsize=(12,5))
plot.show_segmentation(fig, img, masks, flows[0], channels=[0,0])
plt.tight_layout()
plt.show()"""

"""  hyperparameter to play with

    Diameter should we fix it for all image? Changing the diameter will change the results that the algorithm outputs.
    When the diameter is set smaller than the true size then cellpose may over-split cells. Similarly,
    if the diameter is set too big then cellpose may over-merge cells.

    Resample

    3D settings
    flow_threshold: float (optional, default 0.4)
    flow error threshold (all cells with errors below threshold are kept) (not used for 3D)

    cellprob_threshold: float (optional, default 0.0)
    cell probability threshold (all pixels with prob above threshold kept for masks)"""
