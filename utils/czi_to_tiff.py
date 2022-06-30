#!/usr/bin/env python3

#%%

import argparse
import os
from os import listdir
from os.path import isfile, join
import czifile as zis
from matplotlib import pyplot as plt
import tifffile
import numpy as np


def preprare_tiff(path_to_czi, path_to_dapi, path_to_af647, path_to_af568):  # TODO add mkdir
    onlyfiles = [f for f in listdir(path_to_czi) if isfile(join(path_to_czi, f)) and join(path_to_czi, f)[-3:] == "czi"]
    for f in onlyfiles:
        try:
            czi = zis.CziFile(path_to_czi + f)
            metadatadict_czi = czi.metadata(raw=False)     # parse the XML into a dictionary
            chanel = metadatadict_czi['ImageDocument']['Metadata']['DisplaySetting']['Channels']['Channel']
            chanel_name = [chanel[i]["Name"] for i in range(len(chanel))]
            array_im = zis.imread(path_to_czi + f)
            print(array_im.shape)
            for i in range(len(chanel_name)):
                if chanel_name[i] == "Alexa Fluor 647" or chanel_name[i] == 'Cy5':
                    if len(array_im.shape) == 7:
                        array_af647_3d = array_im[0, 0, i, :, :, :, 0]
                    else:
                        array_af647_3d = array_im[0, i, :, :, :, 0]
                    tifffile.imwrite(path_to_af647 + "AF647_" + f[:-3] + "tiff", data=array_af647_3d,
                                     shape=array_af647_3d .shape, dtype=array_af647_3d.dtype)
                if chanel_name[i] == "Alexa Fluor 568" or chanel_name[i] == 'Cy3':
                    if len(array_im.shape) == 7:
                        array_af568_3d = array_im[0, 0, i, :, :, :, 0]
                    else:
                        array_af568_3d = array_im[0, i, :, :, :, 0]
                    tifffile.imwrite(path_to_af568 + "AF568_" + f[:-3] + "tiff", data=array_af568_3d,
                                     shape=array_af568_3d.shape, dtype=array_af568_3d.dtype)
                if chanel_name[i] == "DAPI":
                    if len(array_im.shape) == 7:
                        array_dapi_3d = array_im[0, 0, i, :, :, :, 0]
                    else:
                        array_dapi_3d = array_im[0, i, :, :, :, 0]
                    tifffile.imwrite(path_to_dapi + "dapi_" + f[:-3] + "tiff", data=array_dapi_3d,
                                     shape=array_dapi_3d.shape, dtype=array_dapi_3d.dtype)
        except Exception as e:
            print(e)
            print("ERROR czi")
            print(f)
            print()



def preprare_tiff_for_th_annotation(path_to_czi, path_acquisition = "acquisition"):  # TODO add mkdir
    onlyfiles = [f for f in listdir(path_to_czi) if isfile(join(path_to_czi, f)) and join(path_to_czi, f)[-3:] == "czi"]
    print(onlyfiles)
    print(len(onlyfiles))
    for f in onlyfiles:
        czi = zis.CziFile(path_to_czi + f)
        metadatadict_czi = czi.metadata(raw=False)     # parse the XML into a dictionary
        chanel = metadatadict_czi['ImageDocument']['Metadata']['DisplaySetting']['Channels']['Channel']
        chanel_name = [chanel[i]["Name"] for i in range(len(chanel))]
        array_im = zis.imread(path_to_czi + f)
        t = True
        for i in range(len(chanel_name)):
            image_name = ''.join(e for e in f[:-3] if e.isalnum())
            if chanel_name[i] == "Alexa Fluor 647" or chanel_name[i] == 'Cy5':
                array_af647_3d = array_im[0, i, :, :, :, 0]
                tifffile.imwrite(path_acquisition+ image_name +"_fov0"+"_AF647."  + "tiff", data=array_af647_3d,
                                 shape=array_af647_3d .shape, dtype=array_af647_3d.dtype)
                t = False

            if chanel_name[i] == "Alexa Fluor 568" or chanel_name[i] == 'Cy3':
                array_af568_3d = array_im[0, i, :, :, :, 0]
                tifffile.imwrite(path_acquisition+ image_name + "_fov0"+"_AF568."  + "tiff", data=array_af568_3d,
                                 shape=array_af568_3d.shape, dtype=array_af568_3d.dtype)
                t = False

        if t:
            print(chanel_name)
def f3d_to_mip(path, path_to_save):
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    for f in onlyfiles:
        print(f)
        img = tifffile.imread(path + f)
        img = np.amax(img, 0)
        tifffile.imwrite(path_to_save + f, img)
                               
        
def get_multi_chanel_tiff(path_rna568, path_rna647, path_dapi, path_to_save):
    onlyfiles = [f for f in listdir(path_dapi) if isfile(join(path_dapi, f)) and f[-1] == "f" ]
    onlyfiles = [onlyfiles[i][5:] for i in range(len(onlyfiles))]
    for f in onlyfiles:
        rna_568 = tifffile.imread(path_rna568 + f)
        rna_647 = tifffile.imread(path_rna647 + f)
        nuclei = tifffile.imread(path_dapi +"dapi_" +f) 
        final_array = np.array([rna_568, rna_647, nuclei])
        tifffile.imwrite(path_to_save +  f, final_array)




#%%
if __name__ == "__main__":
    
    
    parser = argparse.ArgumentParser(description='test')
    ### path to datset
    parser.add_argument('-ptz',"--path_to_czi_folder", type=str,
                        default= "/home/tom/Bureau/phd/first_lustra/thesis_sandra_analysis/Lamp3_Pdgfra/NI/NI1225/", 
                        help='path_to_czi folder')

    parser.add_argument('-ptp',"--path_to_project", type=str, 
                        default= "/home/tom/Bureau/phd/first_lustra/thesis_sandra_analysis/Lamp3_Pdgfra/NI/NI1225/",
                        help='path_to_project')
    args = parser.parse_args()


    if not os.path.exists(args.path_to_project + "tiff_data/"):
         os.mkdir(args.path_to_project + "tiff_data/")
    if not os.path.exists(args.path_to_project + "tiff_data/" + "dapi/"):
         os.mkdir(args.path_to_project + "tiff_data/" + "dapi/")
    if not os.path.exists(args.path_to_project + "tiff_data/" + "af568/"):
         os.mkdir(args.path_to_project + "tiff_data/" + "af568/")
    if not os.path.exists(args.path_to_project + "tiff_data/" + "af647/"):
         os.mkdir(args.path_to_project + "tiff_data/" + "af647/")

    path_to_czi = args.path_to_czi_folder
    path_to_dapi = args.path_to_project + "tiff_data/" + "dapi/"
    path_to_af647 = args.path_to_project + "tiff_data/" + "af647/"
    path_to_af568 = args.path_to_project + "tiff_data/" + "af568/"""
    preprare_tiff(path_to_czi, path_to_dapi, path_to_af647, path_to_af568)
    #preprare_tiff_for_th_annotation(path_to_czi = "/home/tom/Bureau/phd/first_lustra/thesis_sandra_analysis/Lamp3_Pdgfra/IRM/IR5M1249/"
    #                               , path_acquisition = "/home/tom/Bureau/phd/first_lustra/thesis_sandra_analysis/Lamp3_Pdgfra_annotation/")



    
    f3d_to_mip("/home/thomas/Bureau/phd/first_lustra/tiff_for_annotation1703/tiff_data/af568/",
     path_to_save = "/home/thomas/Bureau/phd/first_lustra/tiff_for_annotation1703/tiff_data/af568mip/")
    
    
    
    f3d_to_mip("/home/thomas/Bureau/phd/first_lustra/tiff_for_annotation1703/tiff_data/af647/",
     path_to_save = "/home/thomas/Bureau/phd/first_lustra/tiff_for_annotation1703/tiff_data/af647mip/")
    
    f3d_to_mip("//home/thomas/Bureau/phd/first_lustra/tiff_for_annotation1703/tiff_data/dapi/",
     path_to_save = "//home/thomas/Bureau/phd/first_lustra/tiff_for_annotation1703/tiff_data/dapimip/")
    
    
    path_to_af568 = "/home/thomas/Bureau/phd/first_lustra/analyse_image0407/2D_analysis/201030_friyay_charles/tiff_data/af568/"
    path_to_af647 = "/home/thomas/Bureau/phd/first_lustra/analyse_image0407/2D_analysis/201030_friyay_charles/tiff_data/af647/"
    path_acquisition = "/home/thomas/Bureau/phd/first_lustra/analyse_image0407/2D_analysis/201030_friyay_charles/aquisition/"
    
    onlyfiles = [f[6:-5] for f in listdir(path_to_af647) if isfile(join(path_to_af647, f))]
    for f in onlyfiles:
            image_name = ''.join(e for e in f[:] if e.isalnum())
    
            rna647  = tifffile.imread(path_to_af647 + "AF647_" + f + ".tiff")
            rna568  = tifffile.imread(path_to_af568 + "AF568_" + f + ".tiff")
    
            tifffile.imwrite(path_acquisition+ image_name + "_fov0"+"_AF568."  + "tiff", data=rna568,
                                     )
            tifffile.imwrite(path_acquisition+ image_name +"_fov0"+"_AF647."  + "tiff", data=rna647,)