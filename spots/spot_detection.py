#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import time
from os import listdir
from os.path import isfile, join

import bigfish.detection as detection
import bigfish.stack as stack
import numpy as np
import tifffile
from numpy import argmax, nanmax, unravel_index
from scipy.spatial import ConvexHull, Delaunay
from scipy.spatial.distance import pdist, squareform
#%%
from sklearn.cluster import OPTICS, cluster_optics_dbscan


#%%


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

def computer_optics_cluster(spots, eps=25, min_samples = 10, min_cluster_size=10,
                            xi=0.05, scale = np.array([(300/103),1,1]) ):
    """
    Parameters
    apply OPTICS (Ordering Points To Identify the Clustering Structure) https://scikit-learn.org/stable/modules/generated/sklearn.cluster.OPTICS.html
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



