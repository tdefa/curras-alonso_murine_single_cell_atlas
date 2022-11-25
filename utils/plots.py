import colorsys
import os
import sys
import time
from os import listdir
from os.path import isfile, join

import cv2
import czifile as zis
import napari
import networkx as nx
import numpy as np
import tifffile
from matplotlib import pyplot as plt
from matplotlib.patches import RegularPolygon
from scipy import ndimage as ndi
from skimage import data
from skimage.segmentation import find_boundaries

from spots.post_processing import erase_solitary


#%%
def mask_to_rgb(masks, colors = None):
    img = (masks > 0).astype(np.int)
    HSV = np.zeros((img.shape[0], img.shape[1], 3), np.float32)
    HSV[:,:,2] = np.clip(img*1.5, 0, 1.0)
    for n in range(int(masks.max())):
        ipix = (masks==n+1).nonzero()
        if colors is None:
            HSV[ipix[0],ipix[1],0] = np.random.rand()
        else:
            HSV[ipix[0],ipix[1],0] = colors[n,0]
        HSV[ipix[0],ipix[1],1] = 1.0
    RGB = (hsv_to_rgb(HSV) * 255).astype(np.uint8)
    return RGB
def rgb_to_hsv(arr):
    rgb_to_hsv_channels = np.vectorize(colorsys.rgb_to_hsv)
    r, g, b = np.rollaxis(arr, axis=-1)
    h, s, v = rgb_to_hsv_channels(r, g, b)
    hsv = np.stack((h,s,v), axis=-1)
    return hsv

def hsv_to_rgb(arr):
    hsv_to_rgb_channels = np.vectorize(colorsys.hsv_to_rgb)
    h, s, v = np.rollaxis(arr, axis=-1)
    r, g, b = hsv_to_rgb_channels(h, s, v)
    rgb = np.stack((r,g,b), axis=-1)
    return rgb


#%% napari plot

f = "01_NI_CEC-Cy3_Mki67-Cy5_01.tiff"
path_to_project_c = "/media/tom/Elements/to_take/200908_CEC/"
def napari_vizu(path_to_project_c, f, mask = None):
    mask = tifffile.imread(path_to_project_c + "tiff_data/predicted_mask_dapi/dapi_maskdapi_"+f)
    mask = erase_solitary(mask)
    dapi = tifffile.imread(path_to_project_c + "tiff_data/dapi/dapi_" +f)
    af647 = tifffile.imread(path_to_project_c + "tiff_data/af647/AF647_" + f )
    af568 = tifffile.imread(path_to_project_c + "tiff_data/af568/AF568_" + f)
    spots_cy3 = np.load(path_to_project_c + "detected_spot_3d/AF568_"+ f[:-5] +"array.npy")
    spots_cy5 = np.load(path_to_project_c + "detected_spot_3d/AF647_"+ f[:-5] +"array.npy")

    dico_label_cluster = np.load(path_to_project_c + "dico_label_cluster.npy", allow_pickle=True).item()
    [labels_568, labels_647, spots_568, spots_647] = dico_label_cluster[f]
    dico_stat = np.load(path_to_project_c + "dico_stat_2106" + ".npy", allow_pickle=True).item()
    [nb_nuclei, nb_no_rna, nb_cy3, nb_cy5, nb_both, positive_cluster_568, positive_cluster_647, negative_cluster_568, negative_cluster_647] = dico_stat[f]
    nuclei_568_1 = [dico_stat[f][5][i][3] for i in range(len(dico_stat[f][5]))]
    nuclei_647_1 = [dico_stat[f][6][i][3] for i in range(len(dico_stat[f][6]))]
    mask_568 = np.zeros(mask.shape)
    for i in nuclei_568_1:
        mask_568 += mask == i
    mask_647 = np.zeros(mask.shape)
    for i in nuclei_647_1:
        mask_647 += mask == i
    mask_neg = (mask > 0).astype(int) - (mask_647).astype(int) - (mask_568).astype(int)


    cy5_point_class_true = []
    set_cluster_647 = [el[0] for el in positive_cluster_647] + [el[0] for el in positive_cluster_647]
    for s_index in range(len(spots_647)):
        if labels_647[s_index] in set_cluster_647:
            s = spots_647[s_index]
            cy5_point_class_true.append(s)

    cy3_point_class_true = []
    set_cluster_568 = [el[0] for el in positive_cluster_568] + [el[0] for el in negative_cluster_568]
    for s_index in range(len(spots_568)):
        if labels_568[s_index] in set_cluster_568:
            s = spots_568[s_index]
            cy3_point_class_true.append(s)


    viewer = napari.Viewer()
    viewer.add_image(mask_568, name='Pdgfra', scale=(3, 1, 1), rendering = "iso", colormap = "green",iso_threshold = 0.1 )
    viewer.add_image(mask_647, name='lamp3', scale=(3, 1, 1), rendering = "iso", colormap = "red", iso_threshold = 0.1 )
    viewer.add_image(mask_neg, name='mask', scale=(3, 1, 1), rendering = "iso", iso_threshold = 0.1)
    viewer.add_image(dapi, name='dapi', scale=(3, 1, 1), opacity = 0.4)
    viewer.add_image(af568, name='fish_568', scale=(3, 1, 1), opacity= 0.4, colormap = "green")
    viewer.add_image(af647, name='fish_647', scale=(3, 1, 1), opacity = 0.4, colormap = "red")
    viewer.add_points(cy3_point_class_true, name='cy3_point_true', size=3, scale=(3, 1, 1), edge_color="green", face_color="green", ndim=3)
    viewer.add_points(cy5_point_class_true, name='cy5_point_true', size=3, scale=(3, 1, 1), edge_color="red", face_color="red", ndim=3)
    viewer.add_points(spots_568, name='cy3_point', size=3, scale=(3, 1, 1), edge_color="green", face_color="green", symbol = 'ring', ndim=3)
    viewer.add_points(spots_647, name='cy5_point', size=3, scale=(3, 1, 1), edge_color="red", face_color="red", symbol = 'ring', ndim=3)
napari_vizu(path_to_project_c, f, mask = None)