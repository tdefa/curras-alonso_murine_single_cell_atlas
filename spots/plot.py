
#%%
import colorsys
import os
import sys
import time
from os import listdir
from os.path import isfile, join

import alphashape
import cv2
import czifile as zis
import numpy as np
import tifffile
from matplotlib import pyplot as plt
from matplotlib.patches import RegularPolygon
from scipy import ndimage as ndi
from skimage.segmentation import find_boundaries

from spots.post_processing import erase_small_nuclei, erase_solitary

# -*- coding: utf-8 -*-




def mask_image_to_rgb(img, masks, colors = None):

    if img.ndim>2:
        img = img.astype(np.float32).mean(axis=-1)
    else:
        img = img.astype(np.float32)
    img -= img.min()
    img /= img.max()
    HSV = np.zeros((img.shape[0], img.shape[1], 3), np.float32)
    HSV[:,:,2] = np.clip(img*1.5, 0, 1.0)
    for n in range(int(masks.max())):
        ipix = (masks==n+1).nonzero()
        if colors is None:
            HSV[ipix[0],ipix[1],0] = np.random.rand()
        else:
            HSV[ipix[0],ipix[1],0] = colors[n,0]
        HSV[ipix[0],ipix[1],1] = 1.0
    RGB = (hsv_to_rgb(HSV) * 255).astype(np.uint8) #
    return RGB

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


def mask_image_to_rgb2D_from_list(img, masks, nuclei_af568, nuclei_af647, colors = None):
    """

    Args:
        img ():
        masks ():
        nuclei_af568 ():
        nuclei_af647 ():
        colors ():

    Returns:

    """
    if colors is None:
        colors = np.zeros((4,1))
        colors[0,0] = 0.12 #orange
        colors[1,0] = 0.52 #LIGHT BLUE
        colors[2,0] = 0.33 #green
        colors[3,0] = 0.85 #purple
    positive_nuclei  = set(nuclei_af568 + nuclei_af647)
    uncertain_nuclei = set(nuclei_af568) & set(nuclei_af647)
    if img.ndim>2:
        img = np.amax(img, 0).astype(np.float32)
    else:
        img = img.astype(np.float32)
    img -= img.min()
    img /= img.max()
    HSV = np.zeros((img.shape[0], img.shape[1], 3), np.float32)
    HSV[:,:,2] = np.clip(img*1.5, 0, 1.0)
    green = 0
    yellow = 0
    purple = 0
    blue = 0

    for n in np.unique(masks):
        if n==0:
            continue
        ipix = (masks==n).nonzero()
        if n not in positive_nuclei:
            HSV[ipix[0],ipix[1],0] = colors[2,0]
            green += 1

        elif n in uncertain_nuclei:
            HSV[ipix[0],ipix[1],0] = colors[3,0]
            purple  += 1
            HSV[ipix[0],ipix[1],1] = 1.0
        elif n in nuclei_af568:
            HSV[ipix[0],ipix[1],0] = colors[0,0]
            yellow += 1

        else : #it means that n is in nuclei_647
            HSV[ipix[0],ipix[1],0] = colors[1,0]
            blue += 1
        HSV[ipix[0],ipix[1],1] = 1.0


    HSV[:,:,2] = np.clip(img*1.5, 0, 1.0)
    RGB = (hsv_to_rgb(HSV) * 255).astype(np.uint8) #
    return RGB, green, yellow, blue, purple


def mask_image_to_rgb2D_from_list_green_cy3_red_cy5_both_blue_grey(img, masks, nuclei_af568, nuclei_af647, colors = None):
    """

    Args:
        img ():
        masks ():
        nuclei_af568 ():
        nuclei_af647 ():
        colors ():

    Returns:

    """
    if colors is None:
        colors = np.zeros((4,1))
        colors[0,0] = 0.33 #green cy3
        colors[1,0] = 0.01 # red cy5
        colors[2,0] = 0.5 #grey
        colors[3,0] = 0.6 #Blur
    positive_nuclei  = set(nuclei_af568 + nuclei_af647)
    uncertain_nuclei = set(nuclei_af568) & set(nuclei_af647)
    if img.ndim>2:
        img = np.amax(img, 0).astype(np.float32)
    else:
        img = img.astype(np.float32)
    img -= img.min()
    img /= img.max()
    HSV = np.zeros((img.shape[0], img.shape[1], 3), np.float32)
    HSV[:,:,2] = np.clip(img*1.5, 0, 1.0)
    green = 0
    yellow = 0
    purple = 0
    blue = 0

    for n in np.unique(masks):
        if n==0:
            continue
        ipix = (masks==n).nonzero()
        if n not in positive_nuclei:
            HSV[ipix[0],ipix[1],0] = colors[2,0]
            green += 1
            HSV[ipix[0],ipix[1],1] = 0.1
            HSV[ipix[0],ipix[1],2] = np.clip(HSV[ipix[0],ipix[1],2] *1.5, 0, 1)

        elif n in uncertain_nuclei:
            HSV[ipix[0],ipix[1],0] = colors[3,0]
            purple  += 1
            HSV[ipix[0],ipix[1],1] = 0.9
            #HSV[ipix[0],ipix[1],2] = 0.7

        elif n in nuclei_af568:
            HSV[ipix[0],ipix[1],0] = colors[0,0]
            yellow += 1
            HSV[ipix[0],ipix[1],1] = 1

        else : #it means that n is in nuclei_647
            HSV[ipix[0],ipix[1],0] = colors[1,0]
            blue += 1
            HSV[ipix[0],ipix[1],1] = 1


    RGB = (hsv_to_rgb(HSV) * 255).astype(np.uint8) #
    return RGB, green, yellow, blue, purple #green norna, yellow cy3, purle #both  blue #cy5


def mask_image_to_rgb2D_from_list_orange_cy3_other_grey(img, masks, nuclei_af568, nuclei_af647, colors = None):
    if color is None:
        colors = np.zeros((4,1))
        colors[0,0] = 0.12 #orange
        colors[1,0] = 0.01 # red cy5
        colors[2,0] = 0.5 #grey
        colors[3,0] = 0.6 #Blur

    if img.ndim>2:
        img = np.amax(img, 0).astype(np.float32)
    else:
        img = img.astype(np.float32)
    img -= img.min()
    img /= img.max()
    HSV = np.zeros((img.shape[0], img.shape[1], 3), np.float32)
    HSV[:,:,2] = np.clip(img*1.5, 0, 1.0)
    green = 0
    yellow = 0
    purple = 0
    blue = 0

    for n in np.unique(masks):
        if n==0:
            continue
        ipix = (masks==n).nonzero()
        if n not in nuclei_af568:
            HSV[ipix[0],ipix[1],0] = colors[2,0]
            green += 1
            HSV[ipix[0],ipix[1],1] = 0.1
            HSV[ipix[0],ipix[1],2] = np.clip(HSV[ipix[0],ipix[1],2] *1.5, 0, 1)
        else :
            HSV[ipix[0],ipix[1],0] = colors[0,0]
            blue += 1
            HSV[ipix[0],ipix[1],1] = 1
    RGB = (hsv_to_rgb(HSV) * 255).astype(np.uint8) #
    return RGB, green, yellow, blue, purple


def mask_image_to_rgb2D_from_list_orange_cy5_other_grey(img, masks, nuclei_af568, nuclei_af647, colors = None):
    if color is None
        colors = np.zeros((4,1))
        colors[0,0] = 0.12 #orange
        colors[1,0] = 0.01 # red cy5
        colors[2,0] = 0.5 #grey
        colors[3,0] = 0.6 #Blur
    if img.ndim>2:
        img = np.amax(img, 0).astype(np.float32)
    else:
        img = img.astype(np.float32)
    img -= img.min()
    img /= img.max()
    HSV = np.zeros((img.shape[0], img.shape[1], 3), np.float32)
    HSV[:,:,2] = np.clip(img*1.5, 0, 1.0)
    green = 0
    yellow = 0
    purple = 0
    blue = 0

    for n in np.unique(masks):
        if n==0:
            continue
        ipix = (masks==n).nonzero()
        if n not in nuclei_af647:
            HSV[ipix[0],ipix[1],0] = colors[2,0]
            green += 1
            HSV[ipix[0],ipix[1],1] = 0.1
            HSV[ipix[0],ipix[1],2] = np.clip(HSV[ipix[0],ipix[1],2] *1.5, 0, 1)

        else : #it means that n is in nuclei_647
            HSV[ipix[0],ipix[1],0] = colors[0,0]
            blue += 1
            HSV[ipix[0],ipix[1],1] = 1

    RGB = (hsv_to_rgb(HSV) * 255).astype(np.uint8) #
    return RGB, green, yellow, blue, purple


def cluster_over_nuclei_3D(labels, spots, masks, iou_threshold = 0.5, alpha = None):
    positive_cell = []
    nuc_unique = np.sort(np.unique(masks))
    print(np.max(labels))
    for cluster in range(np.max(labels)):
        print(cluster )
        cluster_spots = spots[labels == cluster]
        print(len(cluster_spots))
        if len(cluster_spots) <= 3:
            continue
        t = time.time()
        try:
            t = time.time()
            print("alph")
            alpha_shape = alphashape.alphashape(cluster_spots, alpha)
            print(time.time() - t)
        except Exception as e:
            print(e)
            continue
        for cs in nuc_unique[1:]: #exclude zero
            try:
                t = time.time()

                mask_single_cell = (masks == cs)

                cell_coord = np.array(list(zip(*np.nonzero(mask_single_cell))))
                p1 = np.sum(alpha_shape.contains(cell_coord))

                overlap = p1 / len(cell_coord )
                if overlap > iou_threshold:
                    print((cluster, cs))
                    positive_cell.append(cs)
            except Exception as e:
                    print(e)
    return positive_cell
