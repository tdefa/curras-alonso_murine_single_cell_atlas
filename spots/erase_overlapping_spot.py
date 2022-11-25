# -*- coding: utf-8 -*-

import json
import time

import alphashape
import bigfish.detection as detection
import bigfish.stack as stack
import matplotlib.path as mplPath
import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay, distance
from skimage.draw import polygon

from spots.spot_detection import (cluster_over_nuclei_3D_convex_hull,
                                  computer_optics_cluster)



def erase_overlapping_spot(spots_568, spots_647,
                           kk_568, kk_647,
                           use_z=True, scale=(300, 103, 103)):
    """
    function to erase overlapping detection accros chanel
    Args:
        spots_568 (np.array): spots of the channel Cy3 length wave =568 nm
        spots_647 (np.array):
        kk_568 (int):kk_568 * xy pixel size is the distance radius where co-detection in the Cy3 Cy5 channel are consider as artefact for Cy3,
        kk_647 (int):kk_568 * xy pixel size is the distance radius where co-detection in the Cy3 Cy5 channel are consider as artefact for Cy5
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
    
    if use_z:
        r_dist = np.sqrt(np.square(z_647 * (scale[0]/scale[1]) - z_568.reshape(-1, 1) * (scale[0]/scale[1])) +
                         np.square(x_647 - x_568.reshape(-1, 1)) + np.square(y_647 - y_568.reshape(-1, 1)))
    else:
        r_dist = np.sqrt(np.square(x_647 - x_568.reshape(-1, 1)) + np.square(y_647 - y_568.reshape(-1, 1)))
    new_spots_568 = []
    removed_spots_568 = []
    for s_b in range(len(spots_568)):
        if r_dist[s_b].min() > kk_568:
            new_spots_568.append(spots_568[s_b])
        else:
            removed_spots_568.append(spots_568[s_b])
    new_spots_647 = []
    removed_spots_647 = []
    for s_b in range(len(spots_647)):
        if r_dist[:, s_b].min() > kk_647:
            new_spots_647.append(spots_647[s_b])
        else:
            removed_spots_647.append(spots_647[s_b])
    return np.array(new_spots_568), np.array(removed_spots_568), np.array(new_spots_647), np.array(removed_spots_647)

