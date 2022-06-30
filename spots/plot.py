
#%%
import cv2
import time
import os
from os import listdir
from os.path import isfile, join
import czifile as zis
from matplotlib import pyplot as plt
import tifffile
import numpy as np

from scipy import ndimage as ndi
from skimage.segmentation import find_boundaries

import sys
import colorsys


from matplotlib.patches import RegularPolygon

# -*- coding: utf-8 -*-

from spots.post_processing import erase_solitary, erase_small_nuclei



def mask_image_to_rgb(img, masks, colors = None):
    """if colors is not None:
        if colors.max()>1:
            colors = np.float32(colors)
            colors /= 255
        colors = utils.rgb_to_hsv(colors)"""
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



def plot_nuclei_rna(img, masks, spots, dico_distance, radius = 3,colors = None, framesize=(15, 10), title = 'RNA binding'):
        fig, ax = plt.subplots(1, 2, sharex='col', figsize=framesize)
        fig.suptitle( title, fontsize=16)
        nuclei = mask_image_to_rgb(img, masks, colors = None)
        ax[0].imshow(nuclei)
        ax[1].imshow(nuclei)
        ax[0].set_title('DAPI + nuclei segmentation', fontdict = {'fontsize' : 24})
        ax[1].set_title('nuclei segmentation + RNA detection', fontdict = {'fontsize' : 24})

        for xyz in spots:
            ##define color
            color = nuclei[dico_distance[tuple(xyz)][1]][0] / 255
            x = RegularPolygon((xyz[-1], xyz[-2]), 5, radius, color=color, linewidth=1, fill=True) #or usefct from from matplotlib.patches import RegularPolygon
            ax[1].add_patch(x)
        plt.show()

def plot_dapi_rna(img, masks, spots, dico_distance, radius = 3,colors = None, framesize=(15, 10), title = 'RNA binding'):
        fig, ax = plt.subplots(1, 2, sharex='col', figsize=framesize)
        fig.suptitle( title, fontsize=16)
        nuclei = mask_image_to_rgb(img, masks, colors = None)
        ax[0].imshow(img)
        ax[1].imshow(nuclei)
        ax[0].set_title('DAPI + nuclei segmentation', fontdict = {'fontsize' : 24})
        ax[1].set_title('nuclei segmentation + RNA detection', fontdict = {'fontsize' : 24})

        for xyz in spots:
            ##define color
            color = nuclei[dico_distance[tuple(xyz)][1]][0] / 255
            x = RegularPolygon((xyz[-1], xyz[-2]), 5, radius, color=color, linewidth=1, fill=True) #or usefct from from matplotlib.patches import RegularPolygon
            ax[1].add_patch(x)
        plt.show()




def plot_smfish_rna(img_fish, img_dapi,  masks, spots, dico_distance, radius = 3,colors = None, framesize=(25, 10), title = 'RNA binding'):
        fig, ax = plt.subplots(1, 2, sharex='col', figsize=framesize)
        fig.suptitle( title, fontsize=24)
        nuclei = mask_image_to_rgb(img_dapi, masks, colors = None)

        #MIP of the fish
        img_fish = np.amax(img_fish, 0)

        ax[0].imshow(img_fish)
        ax[1].imshow(img_fish)
        ax[0].set_title('MIP smFISH', fontdict = {'fontsize' : 24})
        ax[1].set_title('MIP smFISH + spot detection', fontdict = {'fontsize' : 24})

        for xyz in spots:
            ##define color
            color = nuclei[dico_distance[tuple(xyz)][1]][0] / 255
            x = RegularPolygon((xyz[-1], xyz[-2]), 5, radius, color=color, linewidth=1, fill=True) #or usefct from from matplotlib.patches import RegularPolygon
            ax[1].add_patch(x)
        plt.show()


def plot_dapi_smfish_rna(img_fish, img_dapi,  masks, spots, dico_distance, radius = 3,colors = None, framesize=(25, 25), title = 'RNA binding'):
        fig, ax = plt.subplots(2, 2, sharex='col', figsize=framesize)
        fig.suptitle( title, fontsize=24)
        nuclei = mask_image_to_rgb(img_dapi, masks, colors = None)

        #MIP of the fish
        img_fish = np.amax(img_fish, 0)

        ax[0, 0].imshow(img_fish)
        ax[0, 1].imshow(img_fish)
        ax[0, 0].set_title('MIP smFISH', fontdict = {'fontsize' : 24})
        ax[0, 1].set_title('MIP smFISH + RNA detection', fontdict = {'fontsize' : 24})
        ax[1, 0].imshow(nuclei)
        ax[1, 1].imshow(nuclei)
        ax[1, 0].set_title('DAPI + nuclei segmentation', fontdict = {'fontsize' : 24})
        ax[1, 1].set_title('nuclei segmentation + RNA detection', fontdict = {'fontsize' : 24})

        for xyz in spots:
            ##define color
            color = nuclei[dico_distance[tuple(xyz)][1]][0] / 255
            x = RegularPolygon((xyz[-1], xyz[-2]), 5, radius, color=color, linewidth=1, fill=True) #or usefct from from matplotlib.patches import RegularPolygon
            ax[0,1].add_patch(x)
            x = RegularPolygon((xyz[-1], xyz[-2]), 5, radius, color=color, linewidth=1, fill=True) #or usefct from from matplotlib.patches import RegularPolygon
            ax[1,1].add_patch(x)

        plt.show()






#%%

if __name__ == "__main__":

    ###plot AF647
    path_to_dapi = "/home/thomas/Bureau/phd/first_one/tiff_data/dapi/"
    path_to_mask_dapi = "/home/thomas/Bureau/phd/first_one/tiff_data/predicted_mask_dapi_st04/"
    path_to_af647 = "/home/thomas/Bureau/phd/first_one/tiff_data/af647/"
    path_to_af568 = "/home/thomas/Bureau/phd/first_one/tiff_data/af568/"


    onlyfiles = [f for f in listdir(path_to_mask_dapi) if isfile(join(path_to_mask_dapi, f)) and f[-1] == "f" ]
    onlyfiles = [onlyfiles[i][14:] for i in range(len(onlyfiles))]

    for f in onlyfiles[:4]:
        img = tifffile.imread(path_to_dapi +"dapi_"+ f)
        img = np.amax(img, 0)

        slice_z = 18
        nuclei = tifffile.imread(path_to_mask_dapi + "dapi_maskdapi_"+ f) # TODO clean path
        erase_small_nuclei(erase_solitary(nuclei), 200)
        img_fish = tifffile.imread(path_to_af647  + "AF647_"+ f) # TODO clean path

        dico_distance = np.load(path_to_af647 + "detected_spot_3d_st04/" + f[:-5] +".npy",allow_pickle='TRUE').item()
        spots = dico_distance.keys()
        #filter spots z
        new_spots = []
        for s in spots:
            if s[0] == slice_z:
                new_spots.append(s)
        plot_nuclei_rna(img, masks = nuclei[slice_z], spots=new_spots, dico_distance=dico_distance, radius = 3,colors = None, framesize=(15, 10), title = f + " AF647")


        #%%
        #plot_smfish_rna(img_fish, img_dapi = img,  masks = nuclei, spots  = spots, dico_distance=dico_distance, radius = 3,colors = None, framesize=(25, 10), title = f)
        plot_nuclei_rna(img, masks = nuclei, spots=spots,
                        dico_distance=dico_distance, radius = 3,colors = None,
                        framesize=(15, 10), title = f + " AF647")
        #plot_dapi_smfish_rna(img_fish,  img_dapi=img,  masks = nuclei, spots=spots, dico_distance=dico_distance, radius = 3,colors = None, framesize=(75, 75), title = f + " AF568")
