# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-


import cv2
import numpy as np
from PIL import Image

import argparse

import cv2
import time
import os
from os import listdir
from os.path import isfile, join
import czifile as zis
from matplotlib import pyplot as plt
import tifffile
import numpy as np
import cellpose
from cellpose import models, io
import bigfish
import bigfish.detection as detection
import bigfish.plot as plot
from scipy import ndimage as ndi
from skimage.segmentation import find_boundaries
from tqdm import tqdm
import random
from post_processing import erase_solitary, erase_small_nuclei
import numpy as np
from pathlib import Path
from skimage.io import imread, imsave
from matplotlib import pyplot as plt
import tifffile
import scipy.stats
from post_processing import erase_solitary, erase_small_nuclei
from skimage.segmentation import mark_boundaries

def evaluate_iou(path_label = "/home/thomas/Bureau/phd/label_dataset/dandra_3d_14",
                 path_prediction = "/home/thomas/Bureau/phd/first_one/tiff_data/predicted_mask_dapi_st04", d3 = True, iou_thresh = 0.5, erase_so = True):
    
    path_label = Path(path_label)
    path_prediction = Path(path_prediction)
    ap_list = []
    ac_list =[]
    for image_files_label in path_label.glob("*/"):
       # print(image_files_label)
        image_label = imread(image_files_label)
        try:
            slice_label = int(list(image_files_label.parts)[-1][-6:-4])
            image_name = str(list(image_files_label.parts)[-1][:-9])
        except :
            slice_label = int(list(image_files_label.parts)[-1][-5:-4])
            image_name = str(list(image_files_label.parts)[-1][:-8])
        if erase_so:
            image_pred =  erase_solitary(tifffile.imread(str(path_prediction) +"/dapi_maskdapi_" +  image_name  + ".tiff"))
        else:
            image_pred =  tifffile.imread(str(path_prediction) +"/dapi_maskdapi_" +  image_name  + ".tiff")
        cell_label = np.unique(image_label)
        cell_pred = np.unique(image_pred[slice_label])
        dico_ap = {}
        for c in cell_label:
            current_cell = (image_label == c).astype(bool)
            for pred in cell_pred:
                overlap = current_cell * (image_pred[slice_label] == pred).astype(bool) # Logical AND
                union = current_cell +  (image_pred[slice_label] == pred).astype(bool)# Logical OR
                if overlap.sum()/float(union.sum()) > iou_thresh:
                    dico_ap[c] = pred
        tp = len(dico_ap)
        fp = len(cell_pred)  - len(dico_ap.values()) #len([cell_pred[i] for i in range(len(cell_pred)) if cell_pred[i] not in dico_ap.values()])
        fn = len(cell_label) - len(dico_ap.keys()) #len([cell_label[i] for i in range(len(cell_label)) if cell_label[i] not in dico_ap.keys()])
    
        ap = tp / (tp+fp+fn)
        ap_list.append(ap)
    return ap_list




if __name__ == "__main__":

    path_label = Path( "/home/thomas/Bureau/phd/label_dataset/dandra_3d_14")
    path_prediction = Path(r"/home/thomas/Bureau/phd/first_one/tiff_data/predicted_mask_dapi_st04")
    path_input = Path(r"/home/thomas/Bureau/phd/first_one/tiff_data/dapi/")
    for image_files_label in path_label.glob("*/"):
        print(image_files_label)
        image_label = imread(image_files_label)
        try:
            slice_label = int(list(image_files_label.parts)[-1][-6:-4])
            image_name = str(list(image_files_label.parts)[-1][:-9])
        except :
            slice_label = int(list(image_files_label.parts)[-1][-5:-4])
            image_name = str(list(image_files_label.parts)[-1][:-8])
        image_pred =  erase_solitary(
        tifffile.imread(str(path_prediction) +"/dapi_maskdapi_" +  image_name  + ".tiff"))
        image_input = tifffile.imread(str(path_input) +"/dapi_" +  image_name  + ".tiff")
        t = mark_boundaries(image_input[slice_label], image_pred[slice_label], color=(0,1,0), mode='large')
        r = mark_boundaries(image_input[slice_label], image_label, color=(1,0,1), mode='large')
        fig, ax = plt.subplots(1,1,  figsize=(130,100))
        ax.imshow(image_input[slice_label],cmap='gray', alpha=1)
        ax.imshow(t, alpha=0.15)
        ax.imshow(r, alpha=0.15)
        plt.show()

