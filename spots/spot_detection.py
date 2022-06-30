#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import argparse

import time
import os
from os import listdir
from os.path import isfile, join
import czifile as zis
from matplotlib import pyplot as plt
import tifffile
import numpy as np
#import cellpose
#from cellpose import models, io

import bigfish
import bigfish.detection as detection
import bigfish.stack as stack
#%%
from spots.post_processing import erase_solitary
from spots.plot import  hsv_to_rgb

from scipy import ndimage as ndi
from utils_ext.utils import get_contours
from spots.post_processing import erase_solitary, erase_small_nuclei
from tqdm import tqdm
#%%
from skimage.transform import resize

from sklearn.cluster import OPTICS, cluster_optics_dbscan
from scipy.spatial.distance import pdist, squareform
from numpy import random, nanmax, argmax, unravel_index



from skimage.segmentation import watershed
import alphashape

from collections import Counter

from scipy.spatial import Delaunay


import matplotlib.path as mplPath
from scipy.spatial import ConvexHull
#%%

### spot

def spot_detection_for_clustering(sigma, rna_path, path_output_segmentaton,
                                  threshold_input = None, 
                                  output_file = "detected_spot_3d/",
                                min_distance = (3,3,3),
                                  path_to_mask_dapi  = None):
    """
    function to detect the spots with a given sigma and return also the theshold
    Parameters
    ----------
    sigma
    float_out
    rna_path
    path_output_segmentaton
    threshold_input
    output_file
    path_to_mask_dapi

    Returns
    -------

    """
    dico_threshold = {}
    onlyfiles = [f for f in listdir(path_output_segmentaton) if isfile(join(path_output_segmentaton, f)) and f[-1] == "f" ]
    onlyfiles = [onlyfiles[i][14:] for i in range(len(onlyfiles))]
    print(onlyfiles)
    for index_path in range(len(rna_path)):
        path = rna_path[index_path]
        for file_index in range(len(onlyfiles)):
            t = time.time()
            rna = tifffile.imread(path + onlyfiles[file_index])
            print(sigma)
            rna_log = stack.log_filter(rna, sigma)#, float_out)
            # local maximum detection
            mask = detection.local_maximum_detection(rna_log, min_distance=min_distance)
            if threshold_input is not None and onlyfiles[file_index] in threshold_input:
                threshold = threshold_input[onlyfiles[file_index]]
                rna_log = stack.log_filter(rna, sigma)#, float_out = False)
                print("manuel threshold")
            else:    
                threshold = detection.automated_threshold_setting(rna_log, mask)
            print(threshold)
            spots, _ = detection.spots_thresholding(rna_log, mask, threshold)
            dico_threshold[onlyfiles[file_index]] = [threshold, len(spots)]
            np.save( output_file + path[-6:] + onlyfiles[file_index][:-5] + 'array.npy', spots)
            print(len(spots))
    return dico_threshold


def spot_detection_for_stiching(sigma, path_to_fish,
                                  threshold_input = None,
                                  output_file = "detected_spot_3d/",
                                min_distance = (3,3,3)):
    """

    Parameters
    ----------
    sigma
    path_to_fish: path to a set of smFish image_ to stich
    threshold_input
    output_file
    min_distance

    Returns
    -------

    """
    dico_threshold = {}
    onlyfiles = [f for f in listdir(path_to_fish) if isfile(join(path_to_fish, f)) and f[-1] == "f" ]
    print(onlyfiles)

    for file_index in range(len(onlyfiles)):
        t = time.time()
        rna = tifffile.imread(path_to_fish + onlyfiles[file_index])
        print(sigma)
        rna_log = stack.log_filter(rna, sigma)#, float_out)
        # local maximum detection
        mask = detection.local_maximum_detection(rna_log, min_distance=min_distance)
        if threshold_input is not None and onlyfiles[file_index] in threshold_input:
            threshold = threshold_input[onlyfiles[file_index]]
            rna_log = stack.log_filter(rna, sigma)#, float_out = False)
            print("manuel threshold")
        else:
            threshold = detection.automated_threshold_setting(rna_log, mask)
        print(threshold)
        spots, _ = detection.spots_thresholding(rna_log, mask, threshold)
        dico_threshold[onlyfiles[file_index]] = [threshold, len(spots)]
        np.save(output_file +  onlyfiles[file_index][:-5] + 'array.npy', spots)
        print(len(spots))
    return dico_threshold

def computer_optics_cluster(spots, eps=25, min_samples = 10, min_cluster_size=10,
                            xi=0.05, scale = np.array([(300/103),1,1]) ):
    """
    Parameters
    ----------
    spots
    eps
    min_samples
    min_cluster_size
    xi
    scale
    Returns label, arrays with the cluster index of each
    -------
    """

    try:
        clust = OPTICS(min_samples=min_samples, xi=xi, min_cluster_size=int(min_cluster_size))
        
        # Run the fit
        if len(scale)==3 and len(spots[0]) ==3:
            print("rescale the clustering")
            print(len(spots))
            clust.fit(spots * scale)
        else:
            clust.fit(spots)        
        labels= cluster_optics_dbscan(reachability=clust.reachability_,
                                           core_distances=clust.core_distances_,
                                               ordering=clust.ordering_, eps=eps)
        return labels
    except ValueError as e:
        print(e)
        return np.array([-1] * len(spots))
    
    
def cluster_over_nuclei_3D_convex_hull(labels, spots, masks, iou_threshold = 0.5, scale = [300, 103, 103], max_spots = 40000):
    """
    Compute the convex hull of each cluster plus cell classification
    Parameters
    ----------
    labels
    spots
    masks
    iou_threshold

    Returns
    -------

    """
    positive_cell = []
    positive_cluster = []
    negative_cluster = []
    mask_single_cell = (masks > 0)
    all_nuclei_coord = np.array(list(zip(*np.nonzero(mask_single_cell)))) #the background is then not taken into account

    if len(spots) > max_spots:
        print(f'abnormaly high number of spots : {len(spots)} ?, return empty list')
        return positive_cell, positive_cluster, negative_cluster
    for cluster in range(np.max(labels)+1):
        print(cluster)
        cluster_spots = spots[labels == cluster]
        print(len(cluster_spots))
        print()
        if len(cluster_spots) <= 3:
            continue
        t = time.time()
        try:
            convex_hull = Delaunay(cluster_spots)
        except Exception as e:
                    print(e)
                    continue
        cluster_spots[:,0] = cluster_spots[:,0] * scale[0] #only to compute the volume point cloud, rescale because z is 3 time xy
        cluster_spots[:,1] = cluster_spots[:,1] * scale[1]#be careful to rescal it only once
        cluster_spots[:,2] = cluster_spots[:,2] * scale[2]
        try:
            D = pdist(cluster_spots)
            D = squareform(D)
        except Exception as e:
            print(e)
            return positive_cell, positive_cluster, negative_cluster
        longuest_distance, [I_row, I_col] = nanmax(D), unravel_index( argmax(D), D.shape )
        all_coord_bool = convex_hull.find_simplex(all_nuclei_coord ) >= 0
        dapi_cordo = all_nuclei_coord.reshape(-1, 3)[all_coord_bool]
        candidate_cells = np.sort(np.unique([masks[tuple(co)] for co in dapi_cordo])) #take only into account cell that ovellap the point cloud
        for cs in candidate_cells: #exclude zero is done before
            try:
                t = time.time()    
                mask_single_cell = (masks == cs)    
                cell_coord = np.array(list(zip(*np.nonzero(mask_single_cell))))
                overlap = np.sum(convex_hull.find_simplex(cell_coord) >= 0) / len(cell_coord) 
                if overlap > iou_threshold: #c'est pas un iou just un threshold
                    #print((cluster, cs))
                    positive_cell.append(cs)
                    print(f"positive cell {cs}")

                    positive_cluster.append([cluster, overlap, ConvexHull(cluster_spots).volume, cs, ConvexHull(cluster_spots).area, longuest_distance, len(cluster_spots)])
                else:
                    if overlap > 0:
                        negative_cluster.append([cluster, overlap, ConvexHull(cluster_spots).volume, cs, len(cluster_spots)])
                        print(f"negative_cluster {cs}")

            except Exception as e:
                    print(e)
        print(time.time()-t)
    return positive_cell, positive_cluster, negative_cluster

def mask_image_to_rgb2D_from_list(img, masks, nuclei_af568, nuclei_af647, colors = None):
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



###cluster function

    
def cluster_in_nuclei(labels, spots, masks, nucleus_threshold = 1):
    positive_cell = []
    nuc_unique = np.unique(masks)
    for cluster in range(np.max(labels)):
        cluster_spots = spots[labels == cluster]
        positive_cluster = []
        for cs in cluster_spots:
            if masks[tuple(cs)] > 0:
                positive_cluster.append(masks[tuple(cs)])
        cluster_res_nuc = dict(Counter(positive_cluster))
        for k in cluster_res_nuc:
            if cluster_res_nuc[k] > nucleus_threshold:
                positive_cell.append(k)
    return positive_cell
                
               
def cluster_over_nuclei2D(labels, spots, masks, iou_threshold = 0.5): 
    positive_cell = []
    nuc_unique = np.unique(masks)
    print(np.max(labels))
    for cluster in range(np.max(labels)):
        print(cluster )
        cluster_spots = spots[labels == cluster]
        print(len(cluster_spots))
        if len(cluster_spots) <= 2:
            continue
        t = time.time()
        grid  = generate_grid(cluster_spots, nx = 1040, ny = 1388)
        print(time.time() - t)
        print()
        masks_cluster = grid * masks
        candidate_cell = np.sort(np.unique(masks_cluster))
        if len(candidate_cell) >= 2:
            for cs in candidate_cell[1:]:
                mask_cs = masks == cs
                masks_cluster_cs = masks_cluster == cs
                overlap = mask_cs * masks_cluster_cs # Logical AND
                union = mask_cs + masks_cluster_cs # Logical OR
                IOU = overlap.sum()/float(union.sum()) 
                if IOU > 0.5:
                    positive_cell.append(cs)
    return positive_cell



def generate_grid(cluster_spots, nx = 1040, ny = 1388):
    try:         
        alpha_shape = alphashape.alphashape(cluster_spots)
        poly_verts = np.array(alpha_shape.exterior.coords).astype(int)
        poly_verts = [[p[1], p[0]] for p in poly_verts]
    
        x, y = np.meshgrid(np.arange(ny), np.arange(nx))
        x, y = x.flatten(), y.flatten()    
        points = np.vstack((x,y)).T
        path = mplPath.Path(poly_verts)
        grid = path.contains_points(points)
        grid = grid.reshape((nx,ny))
        return grid
    except Exception as e:
        print(e)
        return np.zeros([nx, ny])

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
               # print(time.time() - t)
    
                cell_coord = np.array(list(zip(*np.nonzero(mask_single_cell))))
                #print("over")  
                #print(time.time() - t)
                p1 = np.sum(alpha_shape.contains(cell_coord)) 
                #print(time.time() - t)
    
                overlap = p1 / len(cell_coord )
               # print(time.time() - t)
                if overlap > iou_threshold:
                    print((cluster, cs))                    
                    positive_cell.append(cs)
            except Exception as e:
                    print(e)
    return positive_cell




#%%
if __name__ == "__main__":

    """
    #%% kernel approach

    f = "01_NI_Lamp3-Cy5_Pdgfra-Cy3_01.tiff"
    path_to_project_c = "/media/tom/Elements/to_take/200908_fibrosis/"
    dico_label_cluster = np.load(path_to_project_c + "dico_label_cluster.npy", allow_pickle=True).item()
    [labels_568, labels_647, spots_568, spots_647] = dico_label_cluster[f]
    dico_stat = np.load(path_to_project_c + "dico_stat_2106" + ".npy", allow_pickle=True).item()
    from scipy import stats
    kernel = stats.gaussian_kde(spots_647.T)
    coord = []
    for z in range(53):
        for x in range(1048):
            for y in range(1338):
                coord.append([z,x,y])
    def get_density_mask(spots, bandwidith, masks_dim):
    """

    voxel_size_z = 300
    voxel_size_yx = 103
    psf_z = 375
    psf_yx = 129
    parser = argparse.ArgumentParser(description='test')
    ### path to dapi tiff datset
    parser.add_argument("--path_to_mask_dapi", type=str, 
                        
    default= "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Capillary_EC/tiff_data/predicted_mask_dapi/" , help='')
    
    parser.add_argument("--path_to_dapi", type=str,
    default= "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Capillary_EC/tiff_data/dapi/", help='')
   
    parser.add_argument("--path_to_af647", type=str,
    default= "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Capillary_EC/tiff_data/af647/", help='')
    parser.add_argument("--path_to_af568", type=str,
     default= "/home/tom/Bureau/annotation/cell_type_annotation/to_take/200828-NIvsIR5M/00_Capillary_EC/tiff_data/af568/" , help='')
    parser.add_argument("--seg_3d", type=bool, default= False , help='')
    args = parser.parse_args()
    path_to_mask_dapi = args.path_to_mask_dapi
    path_to_af647 = args.path_to_af647
    path_to_af568 = args.path_to_af568
    rna_path = [path_to_af568 + "AF568_", path_to_af647 + "AF647_"]
    onlyfiles = [f for f in listdir(path_to_mask_dapi) if isfile(join(path_to_mask_dapi, f)) and f[-1] == "f" ]
    onlyfiles = [onlyfiles[i][14:] for i in range(len(onlyfiles))]
    
    dapi = tifffile.imread(args.path_to_dapi + "dapi_"  + onlyfiles[0])
    
    mask_dapi = tifffile.imread(path_to_mask_dapi + "dapi_maskdapi_" + onlyfiles[0])
    
    rna = tifffile.imread(rna_path[1] + onlyfiles[0])
    sigma =  (1.25, 1.25, 1.25)
    min_distance = (3,3,3)
    print(sigma)
    # LoG filter
    rna_log = stack.log_filter(rna, sigma) #, float_out = False)

    # local maximum detection
    mask = detection.local_maximum_detection(rna_log, min_distance=min_distance)

    # thresholding
    threshold = detection.automated_threshold_setting(rna_log, mask)
    print("threshold %s" % str(threshold) )
    
    
    spots, _ = detection.spots_thresholding(rna_log, mask, threshold)
    
    fig, ax = plt.subplots(1,1,  figsize=(20,10))
    ax.imshow(np.amax(mask_dapi, 0))

    for s in spots:
        ax.scatter(s[2],s[1], c = "red", s = 5)
    plt.show()

    #%%



    
    #%%

    
#%%
### Old function
    
def rna_nuclei_link(nuclei , spots, voxel_size_z = "e", voxel_size_yx = "e"):
    #if len(nuclei.shape) == 3:
     #nuclei = erase_solitary(nuclei).astype(int)
    dico_distance = {}
    for sp in spots:
        dico_distance[tuple(sp)] = []
    t = time.time()
    dico_nuc = {}  # create a dico with a nucleus and one point inside, useful for plot
    if len(nuclei.shape) == 3:
        pass
        #for slide in tqdm(nuclei):
            #list_nuc_local = sorted(np.unique(slide))[1:]
            #for nuc in list_nuc_local:
                #dico_nuc[nuc] = get_contours(slide, int(nuc), one_point=True) #get_contours seems to work only in 2D ?
    if len(nuclei.shape) == 2:
        list_nuc_local = sorted(np.unique(nuclei))[1:]
        for nuc in list_nuc_local:
            dico_nuc[nuc] = get_contours(nuclei, int(nuc), one_point=True) #get_contours seems to work only in 2D ?

    list_nuc  = sorted(np.unique(nuclei))[1:]
    for nuc in tqdm(list_nuc):
        inverted_mask = np.ones(nuclei.shape) - (nuclei == nuc).astype(np.int)
        if len(nuclei.shape) == 3:
            distance_to_nucleus = ndi.distance_transform_edt(inverted_mask, sampling = [voxel_size_z, voxel_size_yx, voxel_size_yx])
        else:
            distance_to_nucleus = ndi.distance_transform_edt(inverted_mask)  # compute distance map to border
        for sp in spots:
            if len(nuclei.shape) == 3:
                dico_distance[tuple(sp)].append((nuc, distance_to_nucleus[sp[0], sp[1], sp[2]]))
            else:
                dico_distance[tuple(sp)].append((nuc, dico_nuc[nuc], distance_to_nucleus[sp[1], sp[2]]))
    for key in dico_distance.keys():
        dico_distance[key] = min(dico_distance[key], key = lambda t: t[-1])
    return dico_distance


def rna_nuclei_link_watershed(nuclei , spots, voxel_size_z = "e", voxel_size_yx = "e"):
    #if len(nuclei.shape) == 3:
     #nuclei = erase_solitary(nuclei).astype(int)
    dico_distance = {}
    for sp in spots:
        dico_distance[tuple(sp)] = []
    t = time.time()
    dico_nuc = {}  # create a dico with a nucleus and one point inside, useful for plot
    inverted_mask = np.ones(nuclei.shape) - (nuclei != 0).astype(np.int)
    if len(nuclei.shape) == 3:
        distance = ndi.distance_transform_edt(inverted_mask, sampling = [voxel_size_z, voxel_size_yx, voxel_size_yx])
    else:
        distance = ndi.distance_transform_edt(inverted_mask)  # compute distance map to border
    labels = watershed(distance, nuclei)
    if len(nuclei.shape) == 3:
        for sp in spots:
            dico_distance[tuple(sp)] = (labels[sp[0], sp[1], sp[2]], distance[sp[0], sp[1], sp[2]])
    else:
        for sp in spots:
            dico_distance[tuple(sp)] = (labels[sp[1], sp[2]], distance[sp[1], sp[2]])        
    return dico_distance



def count_cell_type(img, masks, dico_repartion_af568, dico_repartion_af647):
    green_neg = []
    yellow_af568 = []
    purple_both = []
    blue_af647 = []
    unique_nuc_list = np.unique(masks.astype(int))
    for n in unique_nuc_list:
        if n==0:
            continue
        if max([dico_repartion_af568[n], dico_repartion_af647[n]]) < 5:
            green.append(n)
        elif dico_repartion_af568[n] / (dico_repartion_af647[n]+ 10e-5) > 2:
            yellow += 1
        elif dico_repartion_af647[n] / (dico_repartion_af568[n]+ 10e-5) > 2:
            blue_af647 += 1
        else:
            purple_both  += 1
    return green_neg, yellow_af568, blue_af647, purple_both



def spot_detection(voxel_size_yx, voxel_size_z, psf_yx, psf_z, path_to_mask_dapi, rna_path,
                   cluster=False, save = True, alpha=0.7, beta=1, output_file = "detected_spot_3d_st04/"):
    onlyfiles = [f for f in listdir(path_to_mask_dapi) if isfile(join(path_to_mask_dapi, f)) and f[-1] == "f" ]
    onlyfiles = [onlyfiles[i][14:] for i in range(len(onlyfiles))]
    for path in rna_path:
        for file_index in range(len(onlyfiles)):
            t = time.time()
            rna = tifffile.imread(path + onlyfiles[file_index])
            spots, threshold = detection.detect_spots(rna, return_threshold=True, voxel_size_z=voxel_size_z,
                                                      voxel_size_yx=voxel_size_yx, psf_z=psf_z, psf_yx=psf_yx)
            if cluster:
                spots, clusters, reference_spot = detection.decompose_cluster(rna, spots, voxel_size_z, voxel_size_yx, psf_z, psf_yx,
                                                                            alpha=alpha,  # alpha impacts the number of spots per clus
                                                                            beta=beta)  # beta impacts the number of detected clusters
            print("detected spots")
            print("\r shape: {0}".format(spots.shape))
            print("\r dtype: {0}".format(spots.dtype))
            print("\r threshold: {0}".format(threshold))
            print(time.time()-t)
            t = time.time()
            # upload the corresponding nucleus mask
            nuclei = tifffile.imread(path_to_mask_dapi + "dapi_maskdapi_"+ onlyfiles[file_index]) # TODO clean path
            #   ###compute distances###

            dico_distance = rna_nuclei_link(nuclei=nuclei, spots = spots, voxel_size_z= voxel_size_z, voxel_size_yx = voxel_size_yx)
            if save:
                np.save(path[:-6] + output_file +onlyfiles[file_index][:-5] + '.npy', dico_distance)
                if cluster:
                    np.save(path[:-6] + output_file +onlyfiles[file_index][:-5]+ "clusters" + '.npy', clusters)
                    np.save(path[:-6] + output_file +onlyfiles[file_index][:-5]+ "spots" + '.npy', spots)
            print(time.time()-t)


def rna_by_cell(dico_spots, nuclei):
    # input dico_spots key rna position, value [nucleus_id, nucleus pos, distance to nucleus]
    # return a dictionary {nuclei_id, number of rna}
    nuclei_list = np.unique(nuclei)
    dico_result = {}
    dico_state = {}
    for nuc in nuclei_list:
        dico_result[nuc] = 0
    for nuc in dico_spots.values():
        dico_result[nuc[0]] += 1
    return dico_result

def mask_image_to_rgb_bis2d(img, masks, dico_repartion_af568, dico_repartion_af647,  colors):
    if img.ndim>2:
         #img = img.astype(np.float32).mean(axis=-1)
        img = np.amax(img, 0).astype(np.float32)
    else:
        img = img.astype(np.float32)
    #img = utils.normalize99(img)
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
        if max([dico_repartion_af568[n], dico_repartion_af647[n]]) < 5:
            HSV[ipix[0],ipix[1],0] = colors[2,0]
            green += 1
        elif dico_repartion_af568[n]/ (dico_repartion_af647[n]+ 10e-5) > 2:
            HSV[ipix[0],ipix[1],0] = colors[0,0]
            yellow += 1
            
        elif dico_repartion_af647[n]/ (dico_repartion_af568[n]+ 10e-5) > 2 :
            HSV[ipix[0],ipix[1],0] = colors[1,0]
            blue += 1
        else:
            HSV[ipix[0],ipix[1],0] = colors[3,0]
            purple  += 1
            HSV[ipix[0],ipix[1],1] = 1.0
    HSV[:,:,2] = np.clip(img*1.5, 0, 1.0)
    RGB = (hsv_to_rgb(HSV) * 255).astype(np.uint8) #
    return RGB, green, yellow, blue, purple


def spot_detection_plusdeep(psf_yx, psf_z, path_to_mask_dapi, rna_path,
                   cluster=False, 
                   save = True, 
                   output_file = "detected_spot_3d_st04/",
                  use_deep = False,
                  list_model = [],
                  erase_solitary_nuc = False):
    onlyfiles = [f for f in listdir(path_to_mask_dapi) if isfile(join(path_to_mask_dapi, f)) and f[-1] == "f" ]
    onlyfiles = [onlyfiles[i][14:] for i in range(len(onlyfiles))]
    print(onlyfiles)
    for index_path in range(len(rna_path)):
        path = rna_path[index_path]
        for file_index in range(len(onlyfiles)):
            t = time.time()
            rna = tifffile.imread(path + onlyfiles[file_index])
            """ spots, threshold = detection.detect_spots(rna,
                                                      return_threshold=True, voxel_size_z=300,
                                                      voxel_size_yx=103, psf_z=psf_z, psf_yx=psf_yx)"""
            sigma =  (1.25, 1.25, 1.25)
            min_distance = (3,3,3)
            print(sigma)
            rna_log = stack.log_filter(rna, sigma)#, float_out=False)
            # local maximum detection
            mask = detection.local_maximum_detection(rna_log, min_distance=min_distance)
            threshold = detection.automated_threshold_setting(rna_log, mask)
            print(threshold)
            spots, _ = detection.spots_thresholding(rna_log, mask, threshold)
            print("spots %s" % len(spots))
            if use_deep:
                print("use deep")
                ##############to complete
                rna_568 = tifffile.imread(rna_path[0] + onlyfiles[0])
                rna_647 = tifffile.imread(rna_path[1] + onlyfiles[0])
                spots = select_real_spot(spots, rna_568, rna_647,
                                         model_artifact = list_model[index_path][0],  offset =16,
                                         transform=list_model[index_path][1],
                                         normalize = list_model[index_path][2])
                print(len(spots))

            nuclei = tifffile.imread(path_to_mask_dapi + "dapi_maskdapi_"+ onlyfiles[file_index]) # TODO clean path
            if erase_solitary_nuc:
                print("number nuc %s" % str(len(np.unique((nuclei)))))
                nuclei = erase_solitary(nuclei)
                print("erase solitary")
                print("number nuc %s" % str(len(np.unique((nuclei)))))


            print(len(spots))
            dico_distance = rna_nuclei_link_watershed(nuclei=nuclei, spots = spots, voxel_size_z= 300, voxel_size_yx = 103)
            np.save( output_file + path[-6:] + onlyfiles[file_index][:-5] + '.npy', dico_distance)



def select_real_spot(spots, rna_568, rna_647, model_artifact,  offset, transform, normalize):
    min_image = min(rna_568.min(),rna_647.min())
    max_image = max(rna_568.max(),rna_647.max())
    new_spots = []
    model_artifact = model_artifact.eval()
    rna_cy3_mip = np.max(rna_568,0)
    rna_cy5_mip = np.max(rna_647,0)
    print(rna_cy5_mip.shape)
    for s in spots:
        t = time.time()
        input_ = np.array([rna_cy3_mip[max(s[1]-offset,0):min(s[1]+offset,1040),
                                       max(s[2]-offset, 0):min(s[2]+offset, 1388)],
                    rna_cy5_mip[max(s[1]-offset,0):min(s[1]+offset,1040),
                                max(s[2]-offset, 0):min(s[2]+offset, 1388)]])
        count_w = 0
        if input_.shape != (2, offset * 2, offset * 2):
            count_w += 1
            #print("array with xrong shape %s" % str(input_.shape))
            #print(s)
            result_temp = np.zeros( (2, offset * 2, offset * 2))
            result_temp[:, :input_.shape[-2],:input_.shape[-1]] = input_
            input_ = result_temp
        input_image = input_
        input_image = resize(input_image , (2, 64 , 64), preserve_range = True)

        #print(count_w)
        if normalize == 1:
            input_image = (input_image - min_image) / (max_image - min_image) #min max normalization
            input_image = torch.tensor(input_image)
            input_image = transform(input_image)
            input_image = torch.unsqueeze(input_image, 0)
        elif normalize == 2:
            input_image = input_image.astype('float64')
            input_image[0] = (input_image[0] - input_image[0].min()) / (input_image[0].max() - input_image[0].min()) #min max normalization
            input_image[1] = (input_image[1] - input_image[1].min()) / (input_image[1].max() - input_image[1].min()) #min max normalization
            input_image[0] = (input_image[0] - input_image[0].mean()) / input_image[0].std() #
            input_image[1] = (input_image[1] - input_image[1].mean()) / input_image[1].std()
            input_image = torch.tensor(input_image)
            input_image = transform(input_image)
            input_image = torch.unsqueeze(input_image, 0)
        outputs = model_artifact(input_image)
        _, preds = torch.max(outputs, 1)
        #print(preds)
        if preds == 0:
            new_spots.append(s)
    return new_spots

def mask_image_to_rgb_bis3d(img, masks, dico_repartion_af568, dico_repartion_af647,  colors):
    if img.ndim>2:
         #img = img.astype(np.float32).mean(axis=-1)
        img = np.amax(img, 0).astype(np.float32)
    else:
        img = img.astype(np.float32)
    #img = utils.normalize99(img)
    img -= img.min()
    img /= img.max()
    if masks.ndim == 3:
        HSV_g = np.zeros(( masks.shape[1], masks.shape[2],3), np.float32)
        HSV_y = np.zeros(( masks.shape[1], masks.shape[2],3), np.float32)
        HSV_p = np.zeros((  masks.shape[1], masks.shape[2],3), np.float32)
        HSV_b = np.zeros(( masks.shape[1], masks.shape[2],3), np.float32)
    green = 0
    yellow = 0
    purple = 0
    blue = 0
    dic_nuc = {}
    dic_nuc['y'] = []
    for n in np.unique(masks):
        if n == 0:
            continue
        ipix = (masks==n).nonzero()
        if max([dico_repartion_af568[n], dico_repartion_af647[n]]) < 5:
            HSV_g[ipix[1], ipix[2], 0 ] = colors[2,0]
            HSV_g[ipix[1], ipix[2],1] = 1.0
    
            green += 1
        elif dico_repartion_af568[n]/ (dico_repartion_af647[n]+ 10e-5) > 2:
            HSV_y[ipix[1], ipix[2], 0]  = colors[0,0]
            HSV_y[ipix[1], ipix[2],1] = 1.0
            dic_nuc['y'].append(n)
            yellow += 1
        elif dico_repartion_af647[n]/ (dico_repartion_af568[n]+ 10e-5) > 2 :
            HSV_b[ipix[1], ipix[2],0]  = colors[1,0]
            HSV_b[ipix[1], ipix[2],1] = 1.0
    
            blue += 1
        else:
            HSV_p[ipix[1], ipix[2], 0] = colors[3,0]
            HSV_p[ipix[1], ipix[2],1] = 1.0
    
            purple  += 1
    
    HSV_g[:,:,2] = 0.7
    RGB_g = (hsv_to_rgb(HSV_g) * 255).astype(np.uint8) #
    HSV_y[:,:,2] = 0.7
    RGB_y = (hsv_to_rgb(HSV_y) * 255).astype(np.uint8) #
    HSV_p[:,:,2] = 0.7
    RGB_p = (hsv_to_rgb(HSV_p) * 255).astype(np.uint8) #
    HSV_b[:,:,2] = 0.7
    RGB_b = (hsv_to_rgb(HSV_b) * 255).astype(np.uint8) #
    print(dic_nuc)
    return [RGB_g, RGB_y, RGB_b, RGB_p ], green, yellow, blue, purple

