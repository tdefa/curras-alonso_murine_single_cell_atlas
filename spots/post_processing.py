# -*- coding: utf-8 -*-


import numpy as np

def erase_solitary(mask): #mask en 3D
    mask_bis = np.zeros(mask.shape)

    ##traiter le cas ou n = 1
    current_nuclei = set(np.unique(mask[0]))
    post_nuclei = set(np.unique(mask[1]))
    nuclei_to_remove =  current_nuclei - post_nuclei
    nuclei_to_keep = current_nuclei - nuclei_to_remove # reminder: set operation are different from arithemtic operation
    for nuc in nuclei_to_keep:
        mask_bis[0] += (mask[0] == nuc) * mask[0]

    for i in range(1, len(mask)-1):
        pre_nuclei = set(np.unique(mask[i-1]))
        current_nuclei = set(np.unique(mask[i]))
        post_nuclei = set(np.unique(mask[i+1]))
        nuclei_to_remove =  current_nuclei - pre_nuclei - post_nuclei
        nuclei_to_keep = current_nuclei - nuclei_to_remove # reminder: set operation are different from arithemtic operation
        for nuc in nuclei_to_keep:
            mask_bis[i] += (mask[i] == nuc) *  mask[i]
    ##traiter le cas ou n = -1
    current_nuclei = set(np.unique(mask[-1]))
    pre_nuclei = set(np.unique(mask[-2]))
    nuclei_to_remove =  current_nuclei - pre_nuclei
    nuclei_to_keep = current_nuclei - nuclei_to_remove # reminder: set operation are different from arithemtic operation
    for nuc in nuclei_to_keep:
        mask_bis[-1] += (mask[-1] == nuc) * mask[-1]
    return mask_bis

def erase_small_nuclei(mask, min_size = 340):
    mask_bis = np.zeros(mask.shape)
    for i in range(len(mask)):
        for nuc in np.unique(mask[i]):
            if np.sum((mask[i] == nuc).astype(int)) > min_size:
                mask_bis[i] += (mask[i] == nuc) * mask[i]
    return mask_bis
