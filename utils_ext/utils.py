# -*- coding: utf-8 -*-


import argparse
import os
import random
import time
from os import listdir
from os.path import isfile, join
from pathlib import Path

import cv2
import czifile as zis
import numpy as np
import tifffile
from matplotlib import pyplot as plt
from PIL import Image
from scipy import ndimage as ndi
from skimage.segmentation import find_boundaries
from tqdm import tqdm


def get_contours(nuclei, nucleus_nb, one_point= True):
    thresh = cv2.inRange(nuclei, nucleus_nb, nucleus_nb)
    result = np.zeros_like(nuclei)
    contours = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    contours = contours[0] if len(contours) == 2 else contours[1]
    cntr = contours[0]
    if one_point is True:
        return cntr[0, :, 1], cntr[0, :, 0] #cntr is inverse
    return cntr




# function from bigfish https://github.com/fish-quant/big-fish/blob/6d3eef93217d2f771d14d329ad9a581b66fbe870/bigfish/detection/spot_detection.py
def _get_candidate_thresholds(pixel_values):
    """Choose the candidate thresholds to test for the spot detection.
    Parameters
    ----------
    pixel_values : np.ndarray
        Pixel intensity values of the image.
    Returns
    -------
    thresholds : np.ndarray, np.float64
        Candidate threshold values.
    """
    # choose appropriate thresholds candidate
    start_range = 0
    end_range = int(np.percentile(pixel_values, 99.9999))
    if end_range < 100:
        thresholds = np.linspace(start_range, end_range, num=100)
    else:
        thresholds = [i for i in range(start_range, end_range + 1)]
    thresholds = np.array(thresholds)

    return thresholds


### soustract noise before laplacian
def soustract_noise(rna_af568, rna_af647):
    rna_af568_n = np.zeros(rna_af568.shape)
    rna_af647_n = np.zeros(rna_af647.shape)
    for c_slice in range(len(rna_af568)):
        rna_af568_n[c_slice] = (rna_af568[c_slice] - rna_af568[c_slice].mean()) / rna_af568[c_slice].std()
        rna_af647_n[c_slice] = (rna_af647[c_slice] - rna_af647[c_slice].mean()) / rna_af647[c_slice].std()

    rna_af568_denoise = rna_af568_n - rna_af647_n
    rna_af568_denoise[rna_af568_denoise <= 0] = 0

    rna_af647_denoise = rna_af647_n - rna_af568_n
    rna_af647_denoise[rna_af647_denoise <= 0] = 0
    return rna_af568_denoise, rna_af647_denoise



#https://stackoverflow.com/questions/902761/saving-a-numpy-array-as-an-image
if __name__ == "__main__":
    list_dico = []
    dico0 ={"path_to_dapi" : "/home/thomas/Bureau/phd/first_lustra/tiff_data/dapi/",
    "path_to_dapi_mip" : "/home/thomas/Bureau/phd/first_one/tiff_data/predicted_mask_dapi_mip/",
    "path_to_save" :"/home/thomas/Bureau/phd/kaibu_data/dapi_rna_0403/",
    "output_name" : "_dapi",
    "d3": False,
    "label":  False,
    "color": 2}
    dico1 = {"path_to_dapi" : "/home/thomas/Bureau/phd/first_lustra/tiff_data/af647/",
    "path_to_save" :"/home/thomas/Bureau/phd/kaibu_data/dapi_rna_0403/",
    "output_name" : "_af647Cy5",
    "d3": False,
    "label":  False,
    "color": 1}
    dico2 = {"path_to_dapi" : "/home/thomas/Bureau/phd/first_lustra/tiff_data/af568/",
    "path_to_save" :"/home/thomas/Bureau/phd/kaibu_data/dapi_rna_0403/",
    "output_name" : "_af568Cy3",
    "d3": False,
    "label":  False,
    "color": 0}
    list_dico = [dico0, dico1, dico2]
    for dico in list_dico:
        onlyfiles = [f for f in listdir(dico['path_to_dapi']) if isfile(join(dico['path_to_dapi'] , f)) and f[-1] == "f" ]
        if  dico['d3'] is True:
                for f in onlyfiles[:]:
                    print(f)
                    img = tifffile.imread(dico['path_to_dapi'] + f)
    
                    slice_ = random.randint(0,53)
                    img = img[slice_]
    
                    tifffile.imwrite(dico['path_to_save'] + f[5:-5] + '_s_' + str(slice_) + dico['output_name'] + ".tiff", img)
                    if label is True:
                        mask = tifffile.imread(dico['path_to_dapi_mip'] + "dapi_maskdapi" + f[4:]).astype(np.uint8)
                        mask = mask[slice_]
                        mask = Image.fromarray(mask)
                        mask.save(path_to_save  +  f[5:-5] + '_s_' + str(slice_) + '_dapi__nuc_label' + ".png")
    
        else:
    
            for f in onlyfiles[:]:
                print(f[6:-5])
                if dico['output_name'] =="_dapi":
                    strart_i = 5
                else:
                    strart_i = 6
    
                img = tifffile.imread(dico['path_to_dapi'] + f,)
                img = np.amax(img, 0)
                img_rgb = np.zeros([img.shape[0], img.shape[1], 3])
                img_rgb[:, :, dico["color"]] = img/img.max()
                tifffile.imwrite(dico['path_to_save'] + f[strart_i :-5] + '_mip' +  dico['output_name']  + ".tiff", img)
                #np.save(dico['path_to_save'] + f[strart_i :-5] + '_mip' +  dico['output_name'], img_rgb)
    
                if dico['label'] is True:
                    mask = tifffile.imread(dico['path_to_dapi_mip'] + "dapi_maskdapi" + f[4:]).astype(np.uint8)
                    mask = Image.fromarray(mask)
                    mask.save(dico['path_to_save']  +  f[strart_i :-5] + '_mip' + '_dapi__nuc_label' + ".png")



#%%
                    
if __name__ == "__main__":
    list_dico = []
    dico0 ={"path_to_dapi" : "/home/thomas/Bureau/phd/first_one/tiff_data/dapi/",
    "path_to_dapi_mip" : "/home/thomas/Bureau/phd/first_one/tiff_data/predicted_mask_dapi_mip/",
    "path_to_save" :"/home/thomas/Bureau/phd/kaibu_data/dapi_dye_rgb_one_image/",
    "prefix" : "dapi_",
    "output_name" : "_dapi",
    "d3": False,
    "label":  True,
    "color": 2}
    dico1 = {"path_to_dapi" : "/home/thomas/Bureau/phd/first_one/tiff_data/af647/",
             "path_to_save" :"/home/thomas/Bureau/phd/kaibu_data/dapi_dye_rgb_one_image/",
    "prefix" : "AF647_",
    "output_name" : "_af647Cy5",
    "d3": False,
    "label":  False,
    "color": 1}
    dico2 = {"path_to_dapi" : "/home/thomas/Bureau/phd/first_one/tiff_data/af568/",
    "path_to_save" :"/home/thomas/Bureau/phd/kaibu_data/dapi_dye_rgb_one_image/",
    "prefix" : "AF568_",
    "output_name" : "_af568Cy3",
    "d3": False,
    "label":  False,
    "color": 0}
    list_dico = [dico0, dico1, dico2]
    dico = dico0
    onlyfiles = [f for f in listdir(dico['path_to_dapi']) if isfile(join(dico['path_to_dapi'] , f)) and f[-1] == "f" ]
    
    
    
    
    for f in onlyfiles[:]:
        img_rgb = np.zeros([1040, 1388, 3])
        for dico in list_dico:
            print(f[6:-5])
            if dico['output_name'] =="_dapi":
                strart_i = 5
            else:
                strart_i = 5
        
            img = tifffile.imread(dico['path_to_dapi'] + dico["prefix" ] + f[strart_i: ])
            img = np.amax(img, 0)
        
            img_rgb[:, :, dico["color"]] = img/img.max()
            #tifffile.imwrite(dico['path_to_save'] + f[strart_i :-5] + '_mip' +  dico['output_name']  + ".tiff", img)
    
            if dico['label'] is True:
                mask = tifffile.imread(dico['path_to_dapi_mip'] + "dapi_maskdapi" + f[4:]).astype(np.uint8)
                mask = Image.fromarray(mask)
                mask.save(dico['path_to_save']  +  f[strart_i :-5] + '_mip' + '_dapi__nuc_label' + ".png")
    
        np.save(dico['path_to_save'] + f[strart_i :-5] + '_mip'+'_cy3_cy5_dapi', img_rgb)
    
