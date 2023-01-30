#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###########################
# Function to run the segmentation with Cellpose
############################

import argparse
from os import listdir
from os.path import isfile, join

import numpy as np
import tifffile
from cellpose import models
from matplotlib import pyplot as plt
from tqdm import tqdm
from utils_ext.cellpose_utilis import stitch3D


def segment_nuclei(path_to_dapi, path_to_mask_dapi, dico_param, model, save=True, ):
    """
    segment dapi image and save  them in th path_to_mask_dapi folder
    Args:
        path_to_dapi (str):
        path_to_mask_dapi (str):
        dico_param (dict):
        model (cellpose modem):
        save (bool):
    Returns:
        None
    """
    onlyfiles = [f for f in listdir(path_to_dapi) if isfile(join(path_to_dapi, f))]
    print(onlyfiles)
    for f in tqdm(onlyfiles):
        print(f)
        img = tifffile.imread(path_to_dapi + f)
        print(img.shape)
        if dico_param["mip"] is True and len(img.shape) == 3:
            img = np.amax(img, 0)
        else:
            if len(img.shape) == 3:
                img = img.reshape(img.shape[0], 1, img.shape[1], img.shape[2])
                print(f'image dapi shape after reshape {img.shape}')
                img = list(img)
        masks, flows, styles, diams = model.eval(img, diameter=dico_param["diameter"],
                                                 channels=[0, 0], flow_threshold=dico_param["flow_threshold"],
                                                 do_3D=dico_param["do_3D"], stitch_threshold=0)
        try:
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
            tifffile.imwrite(path_to_mask_dapi + "dapi_mask" + f , data=masks, dtype=masks.dtype)
            np.save(path_to_mask_dapi + "dico_param.npy", dico_param)
