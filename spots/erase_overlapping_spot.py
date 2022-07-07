# -*- coding: utf-8 -*-


import json
import time

import alphashape
import bigfish.detection as detection
import bigfish.stack as stack
import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay, distance
from skimage.draw import polygon

from spots.spot_detection import (cluster_over_nuclei_3D_convex_hull,
                                  computer_optics_cluster)


def generate_grid(cluster_spots, nx=1040, ny=1388):
    try:
        alpha_shape = alphashape.alphashape(cluster_spots)
        poly_verts = np.array(alpha_shape.exterior.coords).astype(int)
        poly_verts = [[p[1], p[0]] for p in poly_verts]

        x, y = np.meshgrid(np.arange(ny), np.arange(nx))
        x, y = x.flatten(), y.flatten()
        points = np.vstack((x, y)).T
        path = mplPath.Path(poly_verts)
        grid = path.contains_points(points)
        grid = grid.reshape((nx, ny))
        return grid
    except Exception as e:
        print(e)
        return np.zeros([nx, ny])


def erase_overlapping_spot(spots_568, spots_647,
                           kk_568, kk_647, use_z = True, scale = [300, 103, 103]):
    """
    function to erase overlapping detection accros chanel
    Args:
        spots_568 (array): spots of the channel Cy3 length wave =568 nm
        spots_647 (array):
        kk_568 (int): number of voxel radius where co-detection in the Cy3 Cy5 channel are consider as artefact for Cy3
        kk_647 (int):number of voxel radius where co-detection in the Cy3 Cy5 channel are consider as artefact for Cy5
        use_z (bool): True means 3D data
        scale (list):

    Returns:

    """
    print("the scale is ")

    z_568 = np.array([s[0] for s in spots_568])
    x_568 = np.array([s[1] for s in spots_568])
    y_568 = np.array([s[2] for s in spots_568])
    z_647 = np.array([s[0] for s in spots_647])
    x_647 = np.array([s[1] for s in spots_647])
    y_647 = np.array([s[2] for s in spots_647])
    
    if use_z: #to re-check 
        r_dist = np.sqrt(np.square(z_647 * (300/103) - z_568.reshape(-1,1) * (300/103)) +
                         np.square(x_647 - x_568.reshape(-1,1)) + np.square(y_647 - y_568.reshape(-1,1)))
    else:
        r_dist = np.sqrt(np.square(x_647 - x_568.reshape(-1,1)) + np.square(y_647 - y_568.reshape(-1,1)))
    new_spots_568 = []
    removed_spots_568  = []
    for s_b in range(len(spots_568)):
        if r_dist[s_b].min() > kk_568:
            new_spots_568.append(spots_568[s_b])
        else:
            removed_spots_568.append(spots_568[s_b])
    new_spots_647 = []
    removed_spots_647 = []
    for s_b in range(len(spots_647)):
        if r_dist[:,s_b].min() > kk_647:
            new_spots_647.append(spots_647[s_b])
        else:
            removed_spots_647.append(spots_647[s_b])
    return np.array(new_spots_568), np.array(removed_spots_568), np.array(new_spots_647), np.array(removed_spots_647)




def erase_point_in_cluster_2Dalphashape(new_spots, removed_spots, eps=25, min_samples = 5, min_cluster_size=5, xi=0.05,
                                        nx = 1040, ny = 1388):
    """
    compute the alphashape of cluster of erase point from erase_overlapping_spot and erase spot contain in it
    Args:
        new_spots ():
        removed_spots ():
        eps ():
        min_samples ():
        min_cluster_size ():
        xi ():
        nx ():
        ny ():

    Returns:

    """
    removed_spots_2d = np.array([[s[1], s[2]] for s in removed_spots])
    new_spots_2d = np.array([[s[1], s[2]] for s in new_spots])
    try:
        labels = computer_optics_cluster(removed_spots_2d, eps=eps, min_samples = min_samples, 
                                         min_cluster_size=min_cluster_size, xi=xi)
    except Exception as e:
        print(e)
        return new_spots
    if len(labels) == 0:
        print("len label 0")
        return new_spots
    index_to_remove = []
    for cluster in range(np.max(labels)):
            cluster_spots = removed_spots_2d[labels == cluster]
            grid  = generate_grid(cluster_spots, nx = nx, ny = ny)
            for s_index in range(len(new_spots)):
                if grid[tuple(new_spots_2d[s_index])]:
                    index_to_remove.append(s_index)

    new_spots = list(new_spots)
    index_to_remove = list(set(index_to_remove))
    print(np.max(labels))
    print("already removed %s "  % str(len(removed_spots)))
    print("index removed %s "  % str(len(index_to_remove)))
    print()
    for index in sorted(index_to_remove , reverse=True):
        del new_spots[index]
    print("alphashape2D")
    return new_spots

def erase_point_in_cluster_3Dalphashape(new_spots, removed_spots, eps=25, min_samples = 5, min_cluster_size=5, xi=0.05):
    try:
        labels = computer_optics_cluster(removed_spots, eps=eps, min_samples = min_samples, min_cluster_size=min_cluster_size, xi=xi)
    except Exception as e:
        print(e)
        return new_spots
    index_to_remove = []
    
    for cluster in range(np.max(labels)):
        cluster_spots = removed_spots[labels == cluster]
        try:
            t = time.time()
            print("alph")
            alpha_shape = alphashape.alphashape(cluster_spots)
            print(time.time() - t)   
            for s_index in range(len(new_spots)):
                    if alpha_shape.contains([new_spots[s_index]]):
                        index_to_remove.append(s_index)
        except Exception as e:
            print(e)
            continue         
    new_spots = list(new_spots)
    index_to_remove = list(set(index_to_remove))
    print(np.max(labels))
    print("already removed %s "  % str(len(removed_spots)))
    print("index removed %s "  % str(len(index_to_remove)))
    print()
    for index in sorted(index_to_remove , reverse=True):
        del new_spots[index]
    print("alphashape")
    return new_spots




