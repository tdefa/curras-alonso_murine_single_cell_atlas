#!/usr/bin/env python3

#%%

import argparse
import os
from os import listdir
from os.path import isfile, join

import czifile as zis
import numpy as np
import tifffile
from matplotlib import pyplot as plt


def preprare_tiff(path_to_czi, path_to_dapi, path_to_af647, path_to_af568):
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

                               
