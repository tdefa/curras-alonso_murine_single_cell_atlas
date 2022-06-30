#!/usr/bin/env python3

# -*- coding: utf-8 -*-
#%%
import time
import pandas as pd
import argparse
import os
import tifffile
from spots.post_processing import erase_solitary

import time


from scipy import ndimage
###########Function area
#from tvtk.common import configure_input
import numpy as np
#from tvtk.api import tvtk
#from tvtk.common import configure_input
from scipy.spatial.distance import pdist, squareform
from skimage.segmentation import watershed

from scipy import ndimage as ndi#%%
from numpy import random, nanmax, argmax, unravel_index
from skimage.measure import find_contours
from pathlib import Path
import re


#%%

def get_dye(gene_smfish, key_cell_name): #should be move somewhere else
    #print(key_cell_name)
    #print(gene_smfish)
    for gene in gene_smfish:
        if gene+ "-Cy3" in key_cell_name:
            return "Cy3"
        if gene + "-Cy5" in key_cell_name:
            return "Cy5"
        if gene+ "Cy3" in key_cell_name:
            return "Cy3"
        if gene + "Cy5" in key_cell_name:
            return "Cy5"
        if gene in key_cell_name:
            print('the dye is not written =  RISK OF ERROR')
            a = re.search(f'{gene}', key_cell_name)

            remaining_name = key_cell_name[a.end():]

            a = re.search(f'{gene}', key_cell_name)
            remaining_name = key_cell_name[a.end():]
            dico_param_probes = {"Lamp3": (32, 0.42),
                                 "Pdgfra": (35, 0.42),
                                 "Chil3": (20, 0.55),
                                 'Cap': (35, 0.30),
                                 'aCap': (35, 0.30),
                                 'acap': (35, 0.30),  # same but case diff
                                 "Ptprb": (27, 0.45),
                                 "Ptprb1": (27, 0.45),
                                 "Fibin": (27, 0.40),
                                 'C3ar1': (35, 0.45),
                                 'Hhip': (35, 0.25),
                                 'Mki67': (40, 0.30),
                                 "Serpine1": (40, 0.50),
                                 "Apln": (30, 0.40),
                                 "Pecam1": (30, 0.40),
                                 "CEC": (35, 0.30),
                                 "Rtkn2": (27, 0.40),
                                 }
            ### if there are a probes after it means that gene is the first so it is Cy5
            if any(word in remaining_name for word in list(dico_param_probes.keys())):
                return "Cy5"
            #### the second name is  Cy3
            else:
                return "Cy3"
    print(key_cell_name)
    print(gene_smfish)
    raise Exception("Probe not detected")



def compute_contours(binary_mask_3D, scale_z = 300, scale_xy = 103):
    """
    binary_mask_3D: 3D np array binary mask of one nucleus
    scale in nanometre
    compute the boundaries coordianate of the nuclei
    """
    if binary_mask_3D.ndim == 3:
        list_contour =  []
        z, x, y = np.nonzero(binary_mask_3D)
        z_unique = np.unique(z)
        #print(z_unique)
        for z_index in z_unique:
            r = list(find_contours((binary_mask_3D[z_index]).astype(int),level =0.5)[0])
            list_contour += [[z_index, r[i][0], r[i][1]] for i in range(len(r))]
        scale_contours = np.array(list_contour).astype(np.int64)
        scale_contours[:, 0]  = scale_contours[:, 0] * scale_z
        scale_contours[:, 1:]  = scale_contours[:, 1:] * scale_xy
    else : 
        scale_contours = find_contours((binary_mask_3D).astype(int),level =0.5)[0]
    return scale_contours


def compute_shape_index_from_mask_nucleus(binary_mask_3D, scale_z = 300, scale_xy = 103):
    """ compute the shpae index of one nucleus for a binary mask"""
    scale_contours = compute_contours(binary_mask_3D,  scale_z,  scale_xy)
    D = pdist(scale_contours)
    D = squareform(D);
    N, [I_row, I_col] = nanmax(D), unravel_index(argmax(D), D.shape )
    volume = np.sum(binary_mask_3D) * scale_z * scale_xy * scale_xy
    shape_index = volume / ((4/3) *np.pi * ((N/2)**3))
    return shape_index, N, volume


def compute_shape_index_sample(img_dapi_mask, positive_nuclei):
    """
    

    Parameters
    ----------
    img_dapi_mask : TYPE
        DESCRIPTION.
    positive_nuclei : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        shape index of the positive nuclei

    """


    total_shape_index = 0
    for nuc in positive_nuclei:
        shape_index, N, volume = compute_shape_index_from_mask_nucleus((img_dapi_mask==nuc).astype(int), 
                                                 scale_z = np.int64(300), scale_xy = np.int64(103))
        total_shape_index += shape_index
    return total_shape_index / len(positive_nuclei) if len(positive_nuclei) > 0 else 0
    


def compute_average_nuclei_size(img_dapi_mask, positive_nuclei, scale_z = 300, scale_xy = 103, in_micrometre = True):
    total_volume = 0
    for nuc in positive_nuclei:
        total_volume += np.sum((img_dapi_mask == nuc).astype(np.int64)) *  np.int64(scale_z) *  np.int64(scale_xy) *  np.int64(scale_xy)
    if in_micrometre:
        total_volume = total_volume * (10**-9)
    return total_volume / len(positive_nuclei) if len(positive_nuclei) > 0 else "Not define"

def get_experiment_name(sample_name):
    if "NI" in sample_name:
        return sample_name[:5]
    else:
        return sample_name[:7]

def count_positive_cell(dico_stat, key_cell_name, dye, exclude_impossible_solution = True):
    """
    exclude_impossible_solution : parameter to not take into account nuclei that are positive to two incompatible cell type
    dico_stat_value: list with len(np.unique(img_dapi_mask)), nb_no_rna, nb_cy3, nb_cy5, nb_both, 
          positive_cluster_568, positive_cluster_647, negative_cluster_568, negative_cluster_647"""
    cell_cy3 = [c[3] for c in dico_stat[key_cell_name][5]] #all classifye wi
    cell_cy5 = [c[3] for c in dico_stat[key_cell_name][6]]
    print(type(key_cell_name))
    print(key_cell_name)

    
    if exclude_impossible_solution:
        if "Serpine1" in key_cell_name or 'Mki67' in key_cell_name:
            if dye == "Cy3":
                return len( cell_cy3), cell_cy3
            elif dye == "Cy5":
                return len(cell_cy5), cell_cy5
            else:
                raise Exception("Dye not detected")
        elif all(w in key_cell_name for w in ['Pecam1', "Apln"]) or all(w in key_cell_name for w in ['Pecam1', "Ptprb"]) or all(w in key_cell_name for w in ['Hhip', "Pdgfra"]):
            if dye == "Cy3":
                return len(cell_cy3), cell_cy3
            elif dye == "Cy5":
                return len(cell_cy5), cell_cy5
            else:
                raise Exception("Dye not detected")
        else:
            if dye == "Cy3":
                return len([cell for cell in cell_cy3 if cell not in cell_cy5] ), [cell for cell in cell_cy3 if cell not in cell_cy5]                  
            elif dye == "Cy5":
                return len([cell for cell in cell_cy5 if cell not in cell_cy3]), [cell for cell in cell_cy5 if cell not in cell_cy3]
            else:
                raise Exception("Dye not detected")
    else:
        if dye == "Cy3":
                return cell_cy3                
        elif dye == "Cy5":
                return cell_cy5 
        else:
            raise Exception("Dye not detected")
    

def compute_average_size(l_d, nuclei = None):
    """
    l_d : is the list of list (cluster number, overlapp, cluster volume, nuclei)
    nuclei: lsit of int of positive nuclei
    """

    if nuclei is not None and len(nuclei) == 0:
        return "not defined"
    if len(l_d) == 0:
        return "not defined"
    try:
        unique_clusters, frequency_cluster = np.unique(np.array(l_d)[:,0], return_counts=True)
    except Exception as e:
        print(l_d)
        raise e
    dico_int = {}
    if nuclei is not None:
        l_d = [lll for lll in l_d if lll[3] in nuclei]
    for tup in l_d:
        dico_int[tup[0]] = []
    for tup in l_d:
            dico_int[tup[0]].append(tup[2] * (10**(-9))   / frequency_cluster[unique_clusters == tup[0]]   ) #add the point cloud volume of each cluster divide by the number of cell
    total_sum = 0
    for k in dico_int.keys():
        total_sum += np.sum(dico_int[k])
    if len(l_d) > 0:
        return total_sum  / len(l_d)
    return "not defined"   

def compute_average_size_precise(l_d, nuclei = None):
    """
    l_d : is the list of list (cluster number, overlapp, cluster volume, nuclei)
    nuclei: lsit of int of positive nuclei
    take only into account well define point cloud
    """
    if len(l_d) == 0:
        return None, "not defined"   
    
    unique_clusters, frequency_cluster = np.unique(np.array(l_d)[:,0], return_counts=True)

    if nuclei is not None:
        l_d = [lll for lll in l_d if lll[3] in nuclei]
    cluster_freq_1 = unique_clusters[frequency_cluster == 1]

    if len(l_d) == 0:
        return None, "not defined"   
   
    l_d = [lll for lll in l_d if lll[0] in cluster_freq_1 and lll[1] > 0.8]
    
    if len(l_d) == 0:
        return None, "not defined"   
    return np.sum([vol[2]*(10**-9)  for vol in l_d]) /  len(l_d), len(l_d)

def compute_shape_index_area_volume(l_d, nuclei = None):
    """
    l_d : is the list of list (cluster number, overlapp, cluster volume, nuclei, area, longuest_distance)
    nuclei: lsit of int of positive nuclei
    take only into account well define point cloud
    """
    if len(l_d) == 0:
        return None, "not defined"   
    unique_clusters, frequency_cluster = np.unique(np.array(l_d)[:,0], return_counts=True)
    dico_int = {}
    if nuclei is not None:
        l_d = [lll for lll in l_d if lll[3] in nuclei]
    cluster_freq_1 = unique_clusters[frequency_cluster == 1]
    if len(l_d) == 0:
        return None, "not defined"   
    l_d = [lll for lll in l_d if lll[0] in cluster_freq_1 and lll[1] > 0.8]
    shape_index = 0
    for vol_aire in l_d:
        shape_index += vol_aire[2] *3 * np.sqrt(4 * np.pi) / ((vol_aire[4]) * np.sqrt((vol_aire[4])))
    
    if len(l_d) == 0:
        return None, "not defined"   
    return shape_index /  len(l_d), len(l_d)

def compute_shape_index_diag_volume(l_d, nuclei = None):
    """
    

    Parameters
    ----------
    l_d : TYPE
         is the list of list (cluster number, overlapp, cluster volume, nuclei, area, longuest_distance)
    nuclei : TYPE, optional
        DESCRIPTION. Tlsit of int of positive nuclei

    Returns
    -------
    TYPE
            take only into account well define point cloud

    TYPE
        DESCRIPTION.

    """
    if len(l_d) == 0:
        return None, "not defined"
    unique_clusters, frequency_cluster = np.unique(np.array(l_d)[:,0], return_counts=True)
    dico_int = {}
    if nuclei != 0:
        l_d = [lll for lll in l_d if lll[3] in nuclei]
    cluster_freq_1 = unique_clusters[frequency_cluster == 1]
    if len(l_d) == 0:
        return None, "not defined"
    
    l_d = [lll for lll in l_d if lll[0] in cluster_freq_1 and lll[1] > 0.8]
    shape_index = 0
    for vol_aire in l_d:
        shape_index +=  vol_aire[2] / ((4/3) *np.pi * ((vol_aire[5]/2)**3))
    
    if len(l_d) == 0:
        return None, "not defined"
    return shape_index /  len(l_d), len(l_d)
    




def get_k_type_nn(img_dapi_mask, positive_nuclei_source, positive_nuclei_neighbor, 
                  voxel_size_z = 300, voxel_size_yx = 103, labels = None):
    """
    

    Parameters
    ----------
    img_dapi_mask : TYPE
        DESCRIPTION.
    positive_nuclei_source : TYPE
        DESCRIPTION.
    positive_nuclei_neighbor : TYPE
        DESCRIPTION.
    voxel_size_z : TYPE, optional
        DESCRIPTION. The default is 300.
    voxel_size_yx : TYPE, optional
        DESCRIPTION. The default is 103.
    labels : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    score_limit_list : TYPE
        DESCRIPTION.
    score_limit_list_over_neibors : TYPE
        DESCRIPTION.

    """

    score_limit_list = []
    score_limit_list_over_neibors = []
    if labels is None:
        print("labels not found computed on fly")
        inverted_mask = np.ones(img_dapi_mask.shape) - (img_dapi_mask != 0).astype(np.int)
        unique_nuclei_mask = np.unique(img_dapi_mask)
        if len(img_dapi_mask.shape) == 3:
            distance = ndi.distance_transform_edt(inverted_mask, 
                                                  sampling = [300, 103, 103])
        else:
            distance = ndi.distance_transform_edt(inverted_mask)  # compute distance map to border
        t = time.time()
        labels = watershed(distance, img_dapi_mask)
    
    for nucleus_pos in positive_nuclei_source:
        t = time.time()
        tess_curent_nuc = np.max((img_dapi_mask == nucleus_pos).astype(int) * labels)
        #print(tess_curent_nuc)
        frontiers = ndimage.maximum_filter((labels == tess_curent_nuc).astype(int) , size=3)
        neighbors_tess = np.unique(frontiers * labels)[1:] # remove the zero of the background
        #print("neighbors_tess %s" % str(neighbors_tess)) ##tesselation and  segmentation are computed so tessalation numbering is the same than segmentation
        if len(neighbors_tess) > 1:
            score_limit = len([n for n in neighbors_tess if n in set(positive_nuclei_neighbor) and n != nucleus_pos ])
            print(score_limit)
            score_limit_list.append(score_limit)
            score_limit_list_over_neibors.append(score_limit / (len(neighbors_tess) -1)) #minus nucleus_pos, the zero is already removed,  to correct on the cloud 
            print(time.time() -t)
    
        else:
            score_limit_list.append(0)
            score_limit_list_over_neibors.append(0)
            print('no neigbors ?')
            continue
    
    return score_limit_list, score_limit_list_over_neibors


def get_knn_ratio(img_dapi_mask, positive_nuclei_source,
                  positive_nuclei_neighbor):
    """
    

    Parameters
    ----------
    img_dapi_mask : TYPE
        DESCRIPTION.
    positive_nuclei_source : TYPE
        DESCRIPTION.
    positive_nuclei_neighbor : TYPE
        DESCRIPTION.
    Returns
    -------
    knn_ratio_list : List
        DESCRIPTION.
    """
    knn_ratio_list = []
    unique_nuclei_mask = np.unique(img_dapi_mask)
    for nucleus_pos in positive_nuclei_source:
        t = time.time()
        inverted_mask = np.ones(img_dapi_mask.shape) - (img_dapi_mask == nucleus_pos).astype(np.int)
        distance_to_cell = ndi.distance_transform_edt(inverted_mask, sampling = [300, 103, 103])
        sort_neighbors = np.sort([np.min(distance_to_cell[img_dapi_mask == nuc])  for nuc in unique_nuclei_mask if nuc != nucleus_pos]) #list des plus proches voisins
        score_list = []
        for n in sort_neighbors: #quels sont les plus proches voisins du type cherché
            if n in positive_nuclei_neighbor:
                score_list.append(1)
            else:
                score_list.append(0)
        knn_ratio_list.append(knn_ratio_list)
        print(time.time() -t)
    return knn_ratio_list


#%%

def generate_exels_one_cell(list_folder = "",
                            gene_smfish = ['Cap', 'aCap', 'CEC', 'acap'],
                            path_to_take = "/media/tom/Transcend/image_lustra0605/Images_Hugo/",
                            path_save = "/media/tom/Transcend/image_lustra0605/exelsonecell/",
                            dico_stat_name = "dico_seg1005.npy",
                            compute_nuclei_size = False,
                            nuclei_shape_index = False,
                            compute_neighbor = False,
                            scale_z=300, scale_xy=103
                            ):
    """
    Parameters
    ----------
    list_param
    Returns
    -------
    """


    columns= [
            "folder_name",
            "experiment",
            "image_name",
            "gene",
            "dye",
            "nb_nuclei",
            "nb_positive",
            "average_point_cloud_size",
            "average_nuclei_size",
            "average_nuclei_mean_shape_index_volume/diameter",
            "average neighbors of same type",
            "average % of neighbors of same type",
            "average_point cloud " + str(gene_smfish[0]) + " with well defined point cloud only",
            "point cloud "+ str(gene_smfish[0]) + " shape index volume/area",
            "point cloud "+ str(gene_smfish[0]) + " shape index volume/diameter",
            'number of nuclei used for '+ str(gene_smfish[0])+ " estimation",
            ]
    dataframe_per_files_NI = pd.DataFrame(columns=columns)

    dataframe_per_files_IR1M = pd.DataFrame(columns=columns)

    dataframe_per_files_IR2M = pd.DataFrame(columns=columns)

    dataframe_per_files_IR3M = pd.DataFrame(columns=columns)

    dataframe_per_files_IR4M = pd.DataFrame(columns=columns)
    
    dataframe_per_files_IR5M = pd.DataFrame(columns=columns)
    
    dataframe_per_files_other_IRM = pd.DataFrame(columns=columns)

    for folder_name in list_folder:
        """print("the dico_stat_name is hardcoded")
        if folder_name in ["211207_10gy/", "211220_10gy_prolicence/"]:
            dico_stat_name = "dico_100122.npy"
        elif folder_name in [ "210205_Prolicence/gCap senes/",  "210205_Prolicence/aCap prolif/", "210205_Prolicence/aCap senes/"]:
            dico_stat_name = "dico_stat_2810.npy"
        elif folder_name in [ "210205_Prolicence/gCap prolif/"]:
            dico_stat_name = "dico_stat_02022022.npy"
        else:
            dico_stat_name = "dico_stat_2106.npy"""
        print(folder_name)
        dico_stat = np.load(path_to_take + folder_name + dico_stat_name, allow_pickle = True).item()
        onlyfiles = [list(dico_stat.keys())]
        print(onlyfiles)
        sorted_name = np.sort(list(dico_stat.keys()))
        for key_cell_name in sorted_name:
            if not any(word in key_cell_name for word in gene_smfish):
                continue
            t = time.time()
            print(key_cell_name)
            experience_name = get_experiment_name(key_cell_name)
            dye = get_dye(gene_smfish, key_cell_name)
            nb_positive, positive_nuclei = count_positive_cell(dico_stat, key_cell_name, dye)
            average_point_cloud_size = compute_average_size(dico_stat[key_cell_name][5], positive_nuclei) if dye == 'Cy3' else compute_average_size(dico_stat[key_cell_name][6], positive_nuclei) 
            average_nuclei_size = None
            sample_shape_index_nuclei = None
            score_limit_list = []
            score_limit_list_over_neibors = []
            average_precise_point_cloud_size_type = None
            average_shape_index_area_volume_type = None
            average_shape_index_diag_volume_type = None
            nb_c_type = None
            #### shape index computation
            if compute_nuclei_size:
                img_dapi_mask = tifffile.imread(path_output_segmentaton + "dapi_maskdapi_" + key_cell_name)
                img_dapi_mask = erase_solitary(img_dapi_mask)
                average_nuclei_size = compute_average_nuclei_size(img_dapi_mask, positive_nuclei, scale_z=scale_z, scale_xy=scale_xy, in_micrometre=True)
            if nuclei_shape_index:
                print(time.time()-t)
                path_output_segmentaton = path_to_take + folder_name + "tiff_data/" + "predicted_mask_dapi/"

                if not compute_nuclei_size:
                    img_dapi_mask = tifffile.imread(path_output_segmentaton + "dapi_maskdapi_" + key_cell_name)
                    img_dapi_mask = erase_solitary(img_dapi_mask)
                sample_shape_index_nuclei = compute_shape_index_sample(img_dapi_mask, positive_nuclei)
                print(time.time()-t)
                average_shape_index_area_volume_type, nb_c_type = compute_shape_index_area_volume(dico_stat[key_cell_name][5],
                                                                     positive_nuclei) if dye == 'Cy3' else compute_shape_index_area_volume(dico_stat[key_cell_name][6],
                                                                                                            positive_nuclei)
                average_shape_index_diag_volume_type, nb_c_type = compute_shape_index_diag_volume(dico_stat[key_cell_name][5],
                                                                     positive_nuclei) if dye == 'Cy3' else compute_shape_index_diag_volume(dico_stat[key_cell_name][6],
                                                                                                            positive_nuclei)


                average_precise_point_cloud_size_type, nb_c_type = compute_average_size_precise(dico_stat[key_cell_name][5],
                                                                     positive_nuclei) if dye == 'Cy3' else compute_average_size_precise(dico_stat[key_cell_name][6],
                                                                                                            positive_nuclei)

                print("average_shape_index_diag_volume_type %s" % average_shape_index_diag_volume_type)
                print("average_precise_point_cloud_size_type %s" % average_precise_point_cloud_size_type)
                print("nb_c_type %s" % nb_c_type)
            if compute_neighbor:
                labels = None
                score_limit_list, score_limit_list_over_neibors = get_k_type_nn(img_dapi_mask,
                                                                                positive_nuclei_source = positive_nuclei,
                                                        positive_nuclei_neighbor = positive_nuclei,
                                                                                voxel_size_z = 300,
                                                                                voxel_size_yx = 103,
                                                                                labels = labels)



            dico_input_pd = {
                        "folder_name": Path(folder_name).parts[-2:],
                        "experiment":  get_experiment_name(key_cell_name),
                        "image_name":key_cell_name,
                        "gene":gene_smfish,
                        "dye":dye,
                        "nb_nuclei": dico_stat[key_cell_name][0],
                        "nb_positive":nb_positive,
                        "average_point_cloud_size":average_point_cloud_size ,
                        "average_nuclei_size": average_nuclei_size,
                        "average_nuclei_mean_shape_index_volume/diameter": sample_shape_index_nuclei ,
                        "average neighbors of same type": np.mean(score_limit_list) if len(score_limit_list) > 0   else None,
                        "average % of neighbors of same type": np.mean(score_limit_list_over_neibors) if len(score_limit_list_over_neibors) > 0 else None,
                        "average_point cloud " + str(gene_smfish[0]) + " with well defined point cloud only": average_precise_point_cloud_size_type,
                        "point cloud "+ str(gene_smfish[0]) + " shape index volume/area": average_shape_index_area_volume_type,
                        "point cloud "+ str(gene_smfish[0]) + " shape index volume/diameter": average_shape_index_diag_volume_type,
                        'number of nuclei used for '+ str(gene_smfish[0])+ " estimation": nb_c_type,
                        }
            
            if "NI" in experience_name or "Ctrl" in experience_name:
                dataframe_per_files_NI.loc[len(dataframe_per_files_NI)] = dico_input_pd
            elif "IR5M" in experience_name:
                dataframe_per_files_IR5M.loc[len(dataframe_per_files_IR5M)] = dico_input_pd
            elif "IR4M" in experience_name:
                dataframe_per_files_IR4M.loc[len(dataframe_per_files_IR4M)] = dico_input_pd
            elif "IR3M" in experience_name:
                dataframe_per_files_IR3M.loc[len(dataframe_per_files_IR3M)] = dico_input_pd
            elif "IR2M" in experience_name:
                dataframe_per_files_IR2M.loc[len(dataframe_per_files_IR2M)] = dico_input_pd
            elif "IR1M" in experience_name:
                dataframe_per_files_IR1M.loc[len(dataframe_per_files_IR1M)] = dico_input_pd
            else:
                dataframe_per_files_other_IRM.loc[len(dataframe_per_files_other_IRM)] = dico_input_pd


            frames = [dataframe_per_files_NI, 
                     pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
                     dataframe_per_files_other_IRM, 
                     pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
                     dataframe_per_files_IR5M]
            result = pd.concat(frames)
            
            if not os.path.exists(path_save):
                os.mkdir(path_save)
            print(result)
            result.to_pickle(path_save  + str(gene_smfish[0]) +".pkl")
            result.to_excel(path_save  + str(gene_smfish[0]) +'.xls')
            print("save")
            assert len(list(result.columns)) == len(columns)

    frames = [dataframe_per_files_NI,
              pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
              dataframe_per_files_IR1M,
              pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
              dataframe_per_files_IR2M,
              pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
              dataframe_per_files_IR3M,
              pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
              dataframe_per_files_IR4M,
              pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
             dataframe_per_files_other_IRM, 
             pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
             dataframe_per_files_IR5M]
    result = pd.concat(frames)
    result.to_pickle(path_save  + str(gene_smfish[0]) +".pkl")
    result.to_excel(path_save  + str(gene_smfish[0]) +'.xls')
    print("save")
    return result

#%%
def generate_exels_cell_state_type(list_folder,
                                   gene_1,
                                   gene_2,
                                   path_to_take,
                                   path_save,
                                   dico_stat_name,
                                   scale_z = 300,
                                   scale_xy = 103,
                                   compute_morpho_features = False):


    columns = [
            "folder_name",
            "experiment",
            "image_name",
            "gene_cell_type",
            "gene_cell_state",
            
            "nb_nuclei",
            "probe1 only",
            "probe2 only",
            "nb_positive_both",
            
            "average_point_cloud_size (μm3)" + str(gene_1[0]),
            "average well define point_cloud_size (μm3)" + str(gene_1[0]),
            "average_nuclei_size  (μm3)" + str(gene_1[0]),
            "average_nuclei_shape_index" + str(gene_1[0]),
            
            "average_point_cloud_size (μm3)"+ str(gene_2[0]),
            "average well define point_cloud_size (μm3) " + str(gene_2[0]),
            "average_nuclei_size (μm3)"+ str(gene_2[0]),
            "average_nuclei_shape_index "+ str(gene_2[0]),

            "average_point_cloud_size (μm3) positive_to_both",
            "average well define point_cloud_size (μm3) positive_to_both",
            "average_nuclei_size (μm3) positive_to_both",
            "average_nuclei_shape_index_positive_to_both",
            
            "average number of "+ str(gene_2[0])+" neighbors for " + str(gene_1[0]),
            "average % of " + str(gene_2[0]) + " neighbors for" + str(gene_1[0]),
            "average number of "+str(gene_1[0]) +" neighbors for " + str(gene_2[0]),
            "average % of " +str(gene_1[0]) + "neighbors for " +str(gene_2[0]),
            
            "point cloud "+ str(gene_1[0]) + " shape index volume/area",
            "point cloud "+ str(gene_1[0]) + " shape index volume/diameter",
            "point cloud "+ str(gene_2[0]) + " shape index volume/area",
            "point cloud "+ str(gene_2[0]) + " shape index volume/diameter",
            
            "point cloud both shape index volume/area", 
            "point cloud both shape index volume/diameter",
            
            'number of nuclei used for '+ str(gene_1[0]) + "estimation",
            'number of nuclei used for '+ str(gene_2[0]) + "estimation",
             'number of nuclei used for estimation of both ',
             ]
    dataframe_per_files_NI = pd.DataFrame(columns=columns)
    
    dataframe_per_files_IR5M = pd.DataFrame(columns=columns)
    
    dataframe_per_files_other_IRM = pd.DataFrame(columns=columns)

    for folder_name in list_folder:
        """print("the dico_stat_name is hardcoded")
        if folder_name in ["211207_10gy/", "211220_10gy_prolicence/"]:
            dico_stat_name = "dico_100122.npy"
        elif folder_name in [ "210205_Prolicence/gCap senes/",  "210205_Prolicence/aCap prolif/", "210205_Prolicence/aCap senes/"]:
            dico_stat_name = "dico_stat_2810.npy"
        elif folder_name in [ "210205_Prolicence/gCap prolif/"]:
            dico_stat_name = "dico_stat_02022022.npy"
        else:
            dico_stat_name = "dico_stat_2106.npy"
        print(folder_name)"""
        path_output_segmentaton = path_to_take + folder_name + "tiff_data/" + "predicted_mask_dapi/"

    
        dico_stat = np.load(path_to_take + folder_name + dico_stat_name, allow_pickle = True).item()
        onlyfiles = [list(dico_stat.keys())]
            
        sorted_name = np.sort(list(dico_stat.keys()))
        for key_cell_name in sorted_name:
            if not any(word in key_cell_name for word in gene_1):
                continue
            if not any(word in key_cell_name for word in gene_2):
                continue
            t = time.time()
    
            print(key_cell_name)
            local_features = {}

            labels = None
            
            ##########gene_1
            dye_type = get_dye(gene_1, key_cell_name)
            nb_positive_gene1, positive_nuclei_type = count_positive_cell(dico_stat, key_cell_name, dye_type)
            average_point_cloud_size_type = compute_average_size(dico_stat[key_cell_name][5],
                                                                 positive_nuclei_type) if dye_type == 'Cy3' else compute_average_size(
                dico_stat[key_cell_name][6],
                positive_nuclei_type)

            if compute_morpho_features:

                path_output_segmentaton = path_to_take + folder_name + "tiff_data/" + "predicted_mask_dapi/"
                img_dapi_mask = tifffile.imread(path_output_segmentaton + "dapi_maskdapi_" + key_cell_name)
                img_dapi_mask = erase_solitary(img_dapi_mask)
                average_nuclei_size_type = compute_average_nuclei_size(img_dapi_mask, positive_nuclei_type, scale_z = scale_z, scale_xy = scale_xy)
                sample_shape_index_type = compute_shape_index_sample(img_dapi_mask, positive_nuclei_type)
                average_point_cloud_size_type_precise, nb_type = compute_average_size_precise(dico_stat[key_cell_name][5],
                                       positive_nuclei_type) if dye_type == 'Cy3' else compute_average_size_precise(dico_stat[key_cell_name][6], positive_nuclei_type)
                area_vol_shape_index_type, nb_type = compute_shape_index_area_volume(img_dapi_mask, positive_nuclei_type)
                diag_vol_shape_index_type, nb_type = compute_shape_index_diag_volume(img_dapi_mask, positive_nuclei_type)

            
            ###########gene_2
            dye_gene = get_dye(gene_2, key_cell_name)
            nb_positive_gene2, positive_nuclei_state = count_positive_cell(dico_stat, key_cell_name, dye_gene)

            average_point_cloud_size_state = compute_average_size(dico_stat[key_cell_name][5],
                                                                  positive_nuclei_state) if dye_gene == 'Cy3' else compute_average_size(
                dico_stat[key_cell_name][6], positive_nuclei_state)

            if compute_morpho_features:
                average_nuclei_size_state = compute_average_nuclei_size(img_dapi_mask, positive_nuclei_state, scale_z = scale_z, scale_xy = scale_xy)
                sample_shape_index_state = compute_shape_index_sample(img_dapi_mask, positive_nuclei_state)
                average_point_cloud_size_state_precise, nb_state = compute_average_size_precise(dico_stat[key_cell_name][5],
                                      positive_nuclei_state) if dye_type == 'Cy3' else compute_average_size_precise(dico_stat[key_cell_name][6], positive_nuclei_state)
                area_vol_shape_index_state, nb_state = compute_shape_index_area_volume(img_dapi_mask, positive_nuclei_state)
                diag_vol_shape_index_state, nb_state = compute_shape_index_diag_volume(img_dapi_mask, positive_nuclei_state)

            ###########gene_both
        
            positive_both = list(set(positive_nuclei_state) & set( positive_nuclei_type))
            

            
            average_point_cloud_size_both = compute_average_size(dico_stat[key_cell_name][5],
                                                positive_both)

            if compute_morpho_features:
                average_nuclei_size_both = compute_average_nuclei_size(img_dapi_mask, positive_both, scale_z = scale_z, scale_xy = scale_xy)
                sample_shape_index_both = compute_shape_index_sample(img_dapi_mask, positive_both)
                average_point_cloud_size_both_precise, nb_both  = compute_average_size_precise(dico_stat[key_cell_name][5],
                                          positive_both) if dye_type == 'Cy3' else compute_average_size_precise(dico_stat[key_cell_name][6], positive_both)
                area_vol_shape_index_both, nb_both = compute_shape_index_area_volume(img_dapi_mask, positive_both)
                diag_vol_shape_index_both , nb_both = compute_shape_index_diag_volume(img_dapi_mask, positive_both)
            """print(time.time()-t)
            ###############neighbors stat
            try:
                labels = np.load(path_to_labels + "watershelddapi_maskdapi_" + key_cell_name + ".npy")
                print("FILE found for "+ folder_name +"  "+ key_cell_name)

            except:
                print("NO file found for "+ folder_name +"  "+ key_cell_name)
                labels = None
                
            score_limit_list_type_to_state, score_limit_list_over_neibors_type_to_state = get_k_type_nn(img_dapi_mask, positive_nuclei_source = positive_nuclei_type,
                                                                                                        positive_nuclei_neighbor = positive_nuclei_state,  
                                                                                                        voxel_size_z = 300, voxel_size_yx = 103, labels = labels)

            score_limit_list_state_to_type, score_limit_list_over_neibors_state_to_type = get_k_type_nn(img_dapi_mask, positive_nuclei_source =positive_nuclei_state , 
                                                                                                        positive_nuclei_neighbor = positive_nuclei_type, 
                                                                                                        voxel_size_z = 300, voxel_size_yx = 103, labels = labels)

            score_limit_list, score_limit_list_over_neibors = get_k_type_nn(img_dapi_mask,positive_nuclei_source = positive_nuclei, positive_nuclei_neighbor = positive_nuclei, voxel_size_z = 300, voxel_size_yx = 103, labels = labels)
            """
            
            dico_input_pd = {
                
            "folder_name": Path(folder_name).parts[-2:],
            "experiment": get_experiment_name(key_cell_name),
            "image_name": key_cell_name,
            "probe1": gene_1[0],
            "probe2 ": gene_2[0],
            
            "nb_nuclei": dico_stat[key_cell_name][0],
            "probe1 only": nb_positive_gene1,
            "probe2 only": nb_positive_gene2,
            "nb_positive_both": len(positive_both),
            
            "average_point_cloud_size (μm3) " + str(gene_1[0]): average_point_cloud_size_type,
            "average well define point_cloud_size (μm3) " + str(gene_1[0]): None,# average_point_cloud_size_type_precise,
            "average_nuclei_size  (μm3)" + str(gene_1[0]): None, #average_nuclei_size_type ,
            "average_nuclei_shape_index" + str(gene_1[0]): None ,#sample_shape_index_type,
            
            "average_point_cloud_size (μm3)"+ str(gene_2[0]): average_point_cloud_size_state,
            "average well define point_cloud_size (μm3) " + str(gene_2[0]): None ,#average_point_cloud_size_state_precise,
            "average_nuclei_size (μm3)"+ str(gene_2[0]): None, #average_nuclei_size_state,
            "average_nuclei_shape_index "+ str(gene_2[0]): None, # sample_shape_index_state,
            
            "average_point_cloud_size (μm3) positive_to_both": average_point_cloud_size_both, 
            "average well define point_cloud_size (μm3) positive_to_both":None, #average_point_cloud_size_both_precise,
            "average_nuclei_size (μm3) positive_to_both": None, #average_nuclei_size_both,
            "average_nuclei_shape_index_positive_to_both": None, #sample_shape_index_both,
            
            "average number of "+ str(gene_2[0])+" neighbors for " + str(gene_1[0]): None,# np.mean(score_limit_list_type_to_state) if len(score_limit_list_type_to_state) > 0 else None,
            "average % of " + str(gene_2[0]) + " neighbors for" + str(gene_1[0]):  None,# np.mean(score_limit_list_over_neibors_type_to_state) if len(score_limit_list_over_neibors_type_to_state) > 0 else None,
            "average number of "+str(gene_1[0]) +" neighbors for " + str(gene_2[0]): None,# np.mean(score_limit_list_state_to_type) if len(score_limit_list_state_to_type) > 0 else None,
            "average % of " +str(gene_1[0]) + "neighbors for " +str(gene_2[0]): None,# np.mean(score_limit_list_over_neibors_state_to_type) if len(score_limit_list_over_neibors_state_to_type) > 0 else None,
            
            "point cloud "+ str(gene_1[0]) + " shape index volume/area": None,# area_vol_shape_index_type,
            "point cloud "+ str(gene_1[0]) + " shape index volume/diameter": None,# diag_vol_shape_index_type,
            "point cloud "+ str(gene_2[0]) + " shape index volume/area": None,# area_vol_shape_index_state,
            "point cloud "+ str(gene_2[0]) + " shape index volume/diameter": None,# diag_vol_shape_index_state,
            
            "point cloud both shape index volume/area": None,# area_vol_shape_index_both,
            "point cloud both shape index volume/diameter": None,# diag_vol_shape_index_both,
            
            'number of nuclei used for '+ str(gene_1[0]) + "estimation": None,# nb_type,
            'number of nuclei used for '+ str(gene_2[0]) + "estimation": None,# nb_state,
            'number of nuclei used for estimation of both ': None,# nb_both,
            }

             
            experience_name = get_experiment_name(key_cell_name)
            if "NI" in experience_name:
                dataframe_per_files_NI.loc[len(dataframe_per_files_NI)] = dico_input_pd
            elif "IR5M" in experience_name:
                dataframe_per_files_IR5M.loc[len(dataframe_per_files_IR5M)] = dico_input_pd
                
            else:
                dataframe_per_files_other_IRM.loc[len(dataframe_per_files_other_IRM)] = dico_input_pd
                
            frames = [dataframe_per_files_NI, 
            pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
            dataframe_per_files_other_IRM, 
            pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
                     dataframe_per_files_IR5M]
            result = pd.concat(frames)
            if not os.path.exists(path_save):
                os.mkdir(path_save)
            result.to_pickle(path_save  + str(gene_1[0]) +"_"+gene_2[0] +".pkl")
            result.to_excel(path_save  + str(gene_1[0]) + "_" + gene_2[0] +'.xls')


    frames = [dataframe_per_files_NI, 
             pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
             dataframe_per_files_other_IRM, 
             pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
             dataframe_per_files_IR5M]
    result = pd.concat(frames)
    if not os.path.exists(path_save):
        os.mkdir(path_save)
    result.to_pickle(path_save  + str(gene_1[0]) +"_"+gene_2[0] +".pkl")
    result.to_excel(path_save  + str(gene_1[0]) + "_" + gene_2[0] +'.xls')
    print("save")

    return result 







#%%




if __name__ == '__main__':
    #%% one cell
    list_folder = [
        #"211207_10gy/",
        #"211220_10gy_prolicence/",

       # "200709_Round1/",#no data
        #"200710_Round2/" #no data
       # "200828-NIvsIR5M/00_Capillary_EC/", #okd
        #"200828-NIvsIR5M/00_Large_Vessels/", #okd
        #"200828-NIvsIR5M/00_Macrophages/", #okd
        #"200908_CEC/", #okd
        #"200908_fibrosis/", #ok
        #"201030_fridyay/", #ok
        #"201127_AM_fibro/", #okd
        #"210205_Prolicence/aCap_prolif/", #okd
        #"210205_Prolicence/aCap_senes/", #okd
        #"210219_myo_fibros_y_macrophages/",#okd
        #"210412_repeat_fibro/IR5M/", #okd
        #"210412_repeat_fibro/NI/", #okd
        #"210413_rep2/", #okd
        #"210425_angiogenesis/", #ok
        #"210426_repeat3/", #ok
        #"210428_IR5M1236_Lamp3-Cy5_Pdgfra-Cy3/",
     #   "01_lamp3-cy3_normal_01-czi_2021-10-25_1543/",
      #  "02_lamp3-cy3_normal-trueview_02-czi_2021-10-25_1544/",
        #"210205_Prolicence/gCap prolif/", #dico_stat_02022022.npy
        #"210205_Prolicence/gCap senes/", #dico_stat_2810.npy
        #"210205_Prolicence/aCap prolif/",#dico_stat_2810.npy
        #"210205_Prolicence/aCap senes/", #dico_stat_2810.npy
        #"220218_IR3M17gy/",
        #"220317_IR3M10gy_2346_A/",
        #"220318_IR3M10gy_2346_C/",
        #"220325_IR3M17gy_2351_B/",
        #"220425_IR3M10gy_2345_B/",
        #"220502_fibro_gcap_B/",
        #"220304_IR3M17gy_prolicence/",
        #"220317_IR3M10gy_2346_B/",
        #"220324_IR3M17gy_2351_A/",
        #"220422_IR3M10gy_2345_A/",
        #"220429_fibro_gcap_A/",
        #"220503_fibro_gcap_C/",
        #"220505_fibro_gcap_D/",
        #"220506_fibro_gcap_E/",
        "test1/"
    ]

    parser = argparse.ArgumentParser(description='test')


    parser.add_argument("--path_save",
                        type=str,
                        default="/media/tom/Elements1/to_take/test_pipeline/exels/",
                        help='path to save the exels')

    parser.add_argument("--path_folder_of_folders",
                        type=str,
                        default="/media/tom/Elements1/to_take/test_pipeline/",
                        help='path to save the exels')


    parser.add_argument("--list_folder", nargs="+", default=list_folder,  help=' list of folders in the czi folders to analyse') #


    parser.add_argument('--list_probes', type=str, nargs='+', action='append', default=[ ['Lamp3'],  ['Pecam1'],  ['Ptprb'],['Hhip'], ["Rtkn2"],
                 ['Apln'], ['Chil3'],  ['Fibin'], ['C3ar1'],
                  ['Mki67'],  ['Cap', 'aCap', 'CEC', 'acap'],])

    parser.add_argument("--dico_stat_name",
                        type=str,
                        default="finaltest_for_hugo.npy",
                        help='path to save the exels')
    parser.add_argument("--port", default=39949)
    parser.add_argument("--mode", default='client')

    args = parser.parse_args()
    print(args)

    for probes in args.list_probes:
        generate_exels_one_cell(list_folder =args.list_folder,
                            gene_smfish =probes,
                            path_to_take = args.path_folder_of_folders,
                            path_save = args.path_save,
                            dico_stat_name = args.dico_stat_name,
                            compute_nuclei_size = False,
                            nuclei_shape_index = False,
                            compute_neighbor = False)


#%% pair

    """
    list_folder = [
       # "220218_IR3M17gy/",
        #"220317_IR3M10gy_2346_A/",
        #"220318_IR3M10gy_2346_C/",
        #"220325_IR3M17gy_2351_B/",
        #"220425_IR3M10gy_2345_B/",
        #"220502_fibro_gcap_B/",
        #"220304_IR3M17gy_prolicence/",
        #"220317_IR3M10gy_2346_B/",
        #"220324_IR3M17gy_2351_A/",
        #"220422_IR3M10gy_2345_A/",
        #"220429_fibro_gcap_A/",
        "220503_fibro_gcap_C/",
        "220505_fibro_gcap_D/",
        "220506_fibro_gcap_E/",

    ]

    path_to_take = "/media/tom/Transcend/image_lustra0605/Images_200522/"

    list_probes = [ [['Cap', 'aCap', 'CEC',  'acap'], ['Mki67']],
                    [['Cap', 'aCap', 'CEC',  'acap'], ['Serpine1']],

                    [["Fibin"], ['Mki67']],
                    [['Fibin'], ['Serpine1']],

                    [["Apln"], ['Mki67']],
                    [['Apln'], ['Serpine1']],

                    [['Ptprb'], ['Mki67']],

                    [['Ptprb'], ['Serpine1']],

                    [['Pdgfra'], ['Hhip']],

                    [['Pecam1'], ['Apln']],

                    [['Pecam1'], ['Ptprb']],

                    ]

    path_save = "/media/tom/Transcend/image_lustra0605/Images_200522/two_probes/"

    dico_stat = "finaldico_seg2205.npy"
    for probes in list_probes:
        list_param = [list_folder, probes[0],probes[1],  path_to_take, path_save, dico_stat]
        generate_exels_cell_state_type(list_folder = list_folder,
                                       gene_1 = probes[0],
                                       gene_2 = probes[1],
                                       path_to_take = path_to_take,
                                       path_save = path_save,
                                       dico_stat_name = dico_stat,
                                       scale_z=300,
                                       scale_xy=103,
                                       compute_morpho_features=False)



    #%%
    list_folder = [
           "211207_10gy/",
          "211220_10gy_prolicence/",
    ]

    path_to_take = "/home/tom/Bureau/annotation/sandra050122/"

    list_probes = [[['Cap', 'aCap', 'CEC'], ['Mki67']],
                   [['Cap', 'aCap', 'CEC'], ['Serpine1']],

                   [["Fibin"], ['Mki67']],
                   [['Fibin'], ['Serpine1']],

                   [["Apln"], ['Mki67']],
                   [['Apln'], ['Serpine1']],

                   [['Ptprb'], ['Mki67']],
                   [['Ptprb'], ['Serpine1']],

                   [['Pdgfra'], ['Hhip']],

                   [['Pecam1'], ['Apln']],

                   [['Pecam1'], ['Ptprb']],

                   ]

    path_save = "/home/tom/Bureau/annotation/sandra050122/exels_all/pair_cells/"
    dico_stat = "dico_100122.npy"
    for probes in list_probes:
        list_param = [list_folder, probes[0], probes[1], path_to_take, path_save, dico_stat]
        generate_exels_cell_state_type(list_param)
    """