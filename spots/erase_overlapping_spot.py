# -*- coding: utf-8 -*-


import bigfish.stack as stack

import time


from scipy.spatial import distance

import bigfish.detection as detection
#from torchvision import datasets, models, transforms



import alphashape

from scipy.spatial import Delaunay


from matplotlib import pyplot as plt

import json

import numpy as np

from skimage.draw import  polygon

from spots.spot_detection import computer_optics_cluster, generate_grid, cluster_over_nuclei_3D_convex_hull


def erase_overlapping_spot(spots_568, spots_647, kk_568, kk_647, use_z = True): #todo mettre la scale en parametre

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

def geojson_to_label(data_json, img_size, label = "nuclei",  binary_labeling=False):
    """
    Function reading a geojson and returning a label image array

    Args:
      data_json (dict): dictionary in the geojson format containing coordionate of object to label
      img_size (tuple int): size of the original labelled image
      binary_labeling (bool): if True it does not separate connected componnent and use 1 as label for all polygon
      if False the N objects are number 1,2,3,...N
    """

    def make_mask(roi_pos, img_size):
        rr, cc = polygon(roi_pos[:, 0], roi_pos[:, 1])
        # Make sure it's not outside
        rr[rr < 0] = 0
        rr[rr > img_size[0] - 1] = img_size[0] - 1
        cc[cc < 0] = 0
        cc[cc > img_size[1] - 1] = img_size[1] - 1
        # Generate mask
        mask = np.zeros(img_size, dtype=np.uint8)
        mask[rr, cc] = 1
        return mask

    mask_loop = np.zeros(img_size)
    color = 1
    for i in range(len(data_json['features'])):
        if data_json['features'][i]['properties']["label"]== label:
            try:
                reg_pos = np.squeeze(np.asarray(data_json['features'][i]['geometry']['coordinates']))
            except KeyError:
                pass  # GeometryCollection
            if len(reg_pos.shape) == 1:  # case where the label is a point
                reg_pos = np.array([reg_pos])
            reg_pos[:, [0, 1]] = reg_pos[:, [1, 0]]  # inverse coordinate from kaibu
            reg_pos[:, 0] = -1*reg_pos[:, 0]+img_size[0]
            if binary_labeling:
                mask_loop += make_mask(reg_pos, img_size)
            else:
                mask_loop += make_mask(reg_pos, img_size) * (color)  # N objectS are number 1,2,3,...N
            color += 1 
    return mask_loop



    
def compute_rescale_dist(x,y):
    if x.ndim == 2 and y.ndim == 2:
        return distance.euclidean(x, y)
    elif  x.ndim == 3 and y.ndim == 3:
        return distance.euclidean(x, y, [3,1,1])
    else:
        print("ERROR")
        

def compute_pre_rec(dico_pred):
    nb_artifact_total_cy3 =  np.sum([dico_pred[k][0][0] for k in dico_pred.keys()])
    nb_artifact_total_cy5 =  np.sum([dico_pred[k][0][1] for k in dico_pred.keys()])
    nb_real_total_cy3 =  np.sum([dico_pred[k][0][2] for k in dico_pred.keys()])
    nb_real_total_cy5 =  np.sum([dico_pred[k][0][3] for k in dico_pred.keys()])
    
    after_nb_artifact_total_cy3 =  np.sum([dico_pred[k][1][0] for k in dico_pred.keys()])
    after_nb_artifact_total_cy5 =  np.sum([dico_pred[k][1][1] for k in dico_pred.keys()])
    after_nb_real_total_cy3 =  np.sum([dico_pred[k][1][2] for k in dico_pred.keys()])
    after_nb_real_total_cy5 =  np.sum([dico_pred[k][1][3] for k in dico_pred.keys()])
    
    
    tp_cy3 =  np.abs(after_nb_artifact_total_cy3 - nb_artifact_total_cy3)
    p_cy3 = nb_artifact_total_cy3
    fp_cy3 = np.abs(after_nb_real_total_cy3 - nb_real_total_cy3)
    
    n_cy3 = nb_real_total_cy3
    tn_cy3 = after_nb_real_total_cy3 
    fn_cy3 = after_nb_artifact_total_cy3
    
    print(tp_cy3/p_cy3, (tp_cy3 / ( p_cy3 + fp_cy3)))
    recall_cy3 = tp_cy3/p_cy3
    precision_cy3 = tp_cy3 / ( tp_cy3 + fp_cy3)
    
    
    tp_cy5 =  np.abs(after_nb_artifact_total_cy5 - nb_artifact_total_cy5)
    p_cy5 = nb_artifact_total_cy5
    fp_cy5 = np.abs(after_nb_real_total_cy5 - nb_real_total_cy5)
    
    n_cy5= nb_real_total_cy5
    tn_cy5 = after_nb_real_total_cy5 
    fn_cy5 = after_nb_artifact_total_cy5
    
    recall_cy5 = tp_cy5 / p_cy5
    precision_cy5 = tp_cy5 / (tp_cy5 + fp_cy5)    
    return recall_cy3, precision_cy3, recall_cy5, precision_cy5, [p_cy3, tp_cy3, fp_cy3, tn_cy3, fn_cy3], [ p_cy5, tp_cy5, fp_cy5,
                                                                                                             fp_cy5, tn_cy5, fn_cy5]

def compute_pre_rec_mean(dico_pred, mean_bool = False):
    nb_artifact_total_cy3 =  np.array([dico_pred[k][0][0] for k in dico_pred.keys()])
    nb_artifact_total_cy5 =  np.array([dico_pred[k][0][1] for k in dico_pred.keys()])
    nb_real_total_cy3 =  np.array([dico_pred[k][0][2] for k in dico_pred.keys()])
    nb_real_total_cy5 =  np.array([dico_pred[k][0][3] for k in dico_pred.keys()])
    
    after_nb_artifact_total_cy3 =  np.array([dico_pred[k][1][0] for k in dico_pred.keys()])
    after_nb_artifact_total_cy5 =  np.array([dico_pred[k][1][1] for k in dico_pred.keys()])
    after_nb_real_total_cy3 =  np.array([dico_pred[k][1][2] for k in dico_pred.keys()])
    after_nb_real_total_cy5 =  np.array([dico_pred[k][1][3] for k in dico_pred.keys()])
    
    
    tp_cy3 =  np.abs(after_nb_artifact_total_cy3 - nb_artifact_total_cy3)
    p_cy3 = nb_artifact_total_cy3
    fp_cy3 = np.abs(after_nb_real_total_cy3 - nb_real_total_cy3)
    
    n_cy3 = nb_real_total_cy3
    tn_cy3 = after_nb_real_total_cy3
    fn_cy3 = after_nb_artifact_total_cy3
    
    print(tp_cy3/p_cy3, (tp_cy3 / (p_cy3 + fp_cy3)))
    recall_cy3 = tp_cy3/p_cy3
    precision_cy3 = tp_cy3 / ( tp_cy3 + fp_cy3)
    
    
    tp_cy5 =  np.abs(after_nb_artifact_total_cy5 - nb_artifact_total_cy5)
    p_cy5 = nb_artifact_total_cy5
    fp_cy5 = np.abs(after_nb_real_total_cy5 - nb_real_total_cy5)
    
    n_cy5= nb_real_total_cy5
    tn_cy5 = after_nb_real_total_cy5 
    fn_cy5 = after_nb_artifact_total_cy5
    
    recall_cy5 = tp_cy5 / p_cy5
    precision_cy5 = tp_cy5 / (tp_cy5 + fp_cy5)
    
    cy3_stat = [p_cy3, tp_cy3, fp_cy3, n_cy3, tn_cy3, fn_cy3]
    cy5_stat = [p_cy5, tp_cy5, fp_cy5, n_cy5, tn_cy5, fn_cy5]
    
    if mean_bool:
        cy3_stat = [np.mean(el) for el in cy3_stat]
        cy5_stat = [np.mean(el) for el in cy5_stat]

        return np.mean(recall_cy3), np.mean(precision_cy3), np.mean(recall_cy5), np.mean(precision_cy5), cy3_stat, cy5_stat
    else:
        return recall_cy3, precision_cy3, recall_cy5, precision_cy5, cy3_stat, cy5_stat

def compute_pre_rec_mean2(dico_pred):
    for k in dico_pred.keys():
        artifact_cy3  = dico_pred[k][0][0]
        artifact_cy5 = dico_pred[k][0][1]
        real_cy3 =  dico_pred[k][0][2]
        real_total_cy5 =  dico_pred[k][0][3] 
        
        after_artifact_total_cy3 = dico_pred[k][1][0]
        after_artifact_total_cy5 = dico_pred[k][1][1]
        after_real_total_cy3 =  dico_pred[k][1][2] 
        after_real_total_cy5 =  dico_pred[k][1][3] 
        
        
        p_cy3 = len(artifact_cy3)
        p_cy5 = len(artifact_cy5)
 
        set_after_artifact_total_cy3 = set(after_artifact_total_cy3)
        tp_cy3 = len([s for s in artifact_cy3 if s not in set_after_artifact_total_cy3])
        set_after_artifact_total_cy5 = set(after_artifact_total_cy5)
        tp_cy5 = len([s for s in artifact_cy5 if s not in set_after_artifact_total_cy5])      
        
        set_after_real_total_cy3 = set(after_real_total_cy3)
        fp_cy3 = len([s for s in real_cy3 if s not in after_real_total_cy3])
        set_artifact_cy5 = set(artifact_cy5)
        tp_cy5 = len([s for s in artifact_cy5 if s not in after_real_total_cy5])      
        




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

def erase_point_in_cluster_2Dconvex_hull(new_spots, removed_spots, eps=25, min_samples = 5, min_cluster_size=5, xi=0.05):
    removed_spots_2d = np.array([[s[1],s[2]] for s in removed_spots])
    new_spots_2d = np.array([[s[1],s[2]] for s in new_spots])
    try:
        labels = computer_optics_cluster(removed_spots_2d, eps=eps, min_samples = min_samples, min_cluster_size=min_cluster_size, xi=xi)
    except Exception as e:
        print(e)
        return new_spots
    index_to_remove = []
    
    for cluster in range(np.max(labels)):
            cluster_spots = removed_spots_2d[labels == cluster]
            try:
                convex_hull = Delaunay(cluster_spots)
                index_bool  = convex_hull.find_simplex(new_spots_2d)
                ind_true = [i for i, x in enumerate(index_bool) if x]
                index_to_remove += ind_true
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
    return new_spots



def erase_point_in_cluster_3Dconvex_hull(new_spots, removed_spots, eps=25, min_samples = 5, min_cluster_size=5, xi=0.05):
    try:
        labels = computer_optics_cluster(removed_spots, eps=eps, 
                                         min_samples = min_samples, min_cluster_size=min_cluster_size, xi=xi)
    except Exception as e:
        print(e)
        return new_spots
    index_to_remove = []
    
    for cluster in range(np.max(labels)):
            cluster_spots = removed_spots[labels == cluster]
            try:
                convex_hull = Delaunay(cluster_spots)
                index_bool  = convex_hull.find_simplex(new_spots)
                ind_true = [i for i, x in enumerate(index_bool) if x]
                index_to_remove += ind_true
            except Exception as e:
                print(e)
                continue
    print("already removed %s "  % str(len(removed_spots)))
    print("index removed %s "  % str(len(index_to_remove)))
    new_spots = list(new_spots)
    index_to_remove = list(set(index_to_remove))
    print(np.max(labels))
    print("already removed %s "  % str(len(removed_spots)))
    print("index removed %s "  % str(len(index_to_remove)))
    print()
    for index in sorted(index_to_remove , reverse=True):
        del new_spots[index]
    return new_spots
    
#%%


path_to_json3 = ["/home/tom/Bureau/annotation/rna_association/04_IR5M_Chil3-Cy3_Serpine1-Cy5_06_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/03_NI_Chil3-Cy3_Serpine1-Cy5_01/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/03_NI_Chil3-Cy3_Serpine1-Cy5_004_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/04_IR5M_Chil3-Cy3_Serpine1-Cy5_02/target_files_v1/annotation.json",]


onlyfiles3 = [ "04_IR5M_Chil3-Cy3_Serpine1-Cy5_06", "03_NI_Chil3-Cy3_Serpine1-Cy5_01", 
             "03_NI_Chil3-Cy3_Serpine1-Cy5_004",
             "04_IR5M_Chil3-Cy3_Serpine1-Cy5_02",]


#path_to_json3 = [ "/home/tom/Bureau/annotation/rna_association/14_IR5M_Cap-Cy3_Serpine1-Cy5_011.json"]
#onlyfiles3  = [ "14_IR5M_Cap-Cy3_Serpine1-Cy5_011"]
path_to_dapi = "/home/tom/Bureau/annotation/tiff_data2804/dapi/"
path_rna_647 = "/home/tom/Bureau/annotation/tiff_data2804/af647/"
path_rna_568 = "/home/tom/Bureau/annotation/tiff_data2804/af568/"

## state cell
onlyfiles3 = ["12_IR5M_Cap-Cy3_Mki67-Cy5_005",
            "12_IR5M_Cap-Cy3_Mki67-Cy5_008",
            "12_IR5M_Cap-Cy3_Mki67-Cy5_009"]


path_to_json3 = ["/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_005_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_008/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_009_mip/target_files_v2/annotation_c.json"
]

onlyfiles3 = [ "08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_010",
 "07_CtrlNI_Pdgfra-Cy3_Serpine1-Cy5_002"]

path_to_json3 = [ "/home/tom/Bureau/annotation/rna_association/08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_010_mip/target_files_v2/annotation.json",
"/home/tom/Bureau/annotation/rna_association/07_CtrlNI_Pdgfra-Cy3_Serpine1-Cy5_002.json"]





onlyfiles3 = ["11_NI_Cap-Cy3_Mki67-Cy5_002",
              "12_IR5M_Cap-Cy3_Mki67-Cy5_005",
            "12_IR5M_Cap-Cy3_Mki67-Cy5_008",
            "12_IR5M_Cap-Cy3_Mki67-Cy5_009"]


path_to_json3 = [ "/home/tom/Bureau/annotation/rna_association/11_NI_Cap-Cy3_Mki67-Cy5_002_mip/target_files_v1/annotation.json",

    "/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_005_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_008/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_009_mip/target_files_v2/annotation_c.json"
]

path_to_json3 = ["/home/tom/Bureau/annotation/rna_association/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06.json", 
"/home/tom/Bureau/annotation/rna_association/01_NI_Lamp3-Cy5_Pdgfra-Cy3_01_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10_mip/target_files_v3/annotation.json"]


onlyfiles3  = ["03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_06", "01_NI_Lamp3-Cy5_Pdgfra-Cy3_01",  "03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10"]



path_to_json3 = ["/home/tom/Bureau/annotation/rna_association/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06.json", 
"/home/tom/Bureau/annotation/rna_association/01_NI_Lamp3-Cy5_Pdgfra-Cy3_01_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10_mip/target_files_v3/annotation.json"]


onlyfiles3  = ["03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_06", "01_NI_Lamp3-Cy5_Pdgfra-Cy3_01",  "03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10"]


onlyfiles1 = ["01_NI_Chil3-Cy3_Mki67-Cy5_02",
             "04_IR5M_Chil3-Cy3_Serpine1-Cy5_06", 
             "08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_010",
             "11_NI_Cap-Cy3_Mki67-Cy5_002",
             "07_CtrlNI_Pdgfra-Cy3_Serpine1-Cy5_002", 
             "14_IR5M_Cap-Cy3_Serpine1-Cy5_011",
             "03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_06",
              "01_NI_Lamp3-Cy5_Pdgfra-Cy3_01",
             "03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10"]

path_to_json1 = [
    "/home/tom/Bureau/annotation/rna_association/01_NI_Chil3-Cy3_Mki67-Cy5_02_mip/target_files_v1/annotation.json", 
"/home/tom/Bureau/annotation/rna_association/04_IR5M_Chil3-Cy3_Serpine1-Cy5_06_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_010_mip/target_files_v2/annotation.json",
 "/home/tom/Bureau/annotation/rna_association/11_NI_Cap-Cy3_Mki67-Cy5_002_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/07_CtrlNI_Pdgfra-Cy3_Serpine1-Cy5_002.json",
"/home/tom/Bureau/annotation/rna_association/14_IR5M_Cap-Cy3_Serpine1-Cy5_011.json",
 "/home/tom/Bureau/annotation/rna_association/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06.json", 
"/home/tom/Bureau/annotation/rna_association/01_NI_Lamp3-Cy5_Pdgfra-Cy3_01_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10_mip/target_files_v3/annotation.json"]


onlyfiles2 = ["02_IR5M_Chil3-Cy3_Mki67-Cy5_05",
             "03_NI_Chil3-Cy3_Serpine1-Cy5_01", 
             "03_NI_Chil3-Cy3_Serpine1-Cy5_004",
             "04_IR5M_Chil3-Cy3_Serpine1-Cy5_02",
            "12_IR5M_Cap-Cy3_Mki67-Cy5_005",
            "12_IR5M_Cap-Cy3_Mki67-Cy5_008",
            "12_IR5M_Cap-Cy3_Mki67-Cy5_009"
    ]

path_to_json2 = [
"/home/tom/Bureau/annotation/rna_association/02_IR5M_Chil3-Cy3_Mki67-Cy5_05_mip/target_files_v1/annotation.json", 
"/home/tom/Bureau/annotation/rna_association/03_NI_Chil3-Cy3_Serpine1-Cy5_01/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/03_NI_Chil3-Cy3_Serpine1-Cy5_004_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/04_IR5M_Chil3-Cy3_Serpine1-Cy5_02/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_005_mip/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_008/target_files_v1/annotation.json",
"/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_009_mip/target_files_v2/annotation_c.json"
]


path_to_json3 = path_to_json1 + path_to_json2 
onlyfiles3  = onlyfiles1 + onlyfiles2



    #%%% spots detection if needed
    
if __name__ == "__main__":

    
    ########af568
    
    for f in onlyfiles3:
        rna = tifffile.imread(path_rna_568+ "AF568_"+ f + '.tiff')
        sigma =  (1.25, 1.25, 1.25)
        min_distance = (3, 3, 3)
        print(sigma)
        rna_log = stack.log_filter(rna, sigma, float_out=False)
        mask = detection.local_maximum_detection(rna_log, min_distance=min_distance)
        threshold = detection.automated_threshold_setting(rna_log, mask)
        print(threshold)

        spots, _ = detection.spots_thresholding(rna_log, mask, threshold)
        print(len(spots))
        np.save("/home/tom/Bureau/annotation/tiff_data2804/spots_af568/" +f, spots)
    for f in onlyfiles3:
        rna = tifffile.imread(path_rna_647+ "AF647_"+ f + '.tiff')
        sigma =  (1.35, 1.35, 1.35)
        min_distance = (3, 3, 3)
        print(sigma)
  

            
        rna_log = stack.log_filter(rna, sigma, float_out=True)
        # local maximum detection
        mask = detection.local_maximum_detection(rna_log, min_distance=min_distance)
        threshold = detection.automated_threshold_setting(rna_log, mask)
        print(threshold)
        spots, _ = detection.spots_thresholding(rna_log, mask, threshold)
        if f == '04_IR5M_Chil3-Cy3_Serpine1-Cy5_06':
            rna_log = stack.log_filter(rna, sigma, float_out=False)
            mask = detection.local_maximum_detection(rna_log, min_distance=min_distance)

            spots, _ = detection.spots_thresholding(rna_log, mask, 4)
            
        np.save("/home/tom/Bureau/annotation/tiff_data2804/spots_af647/" +f, spots)
        
        
#%%
if __name__ == "__main__":
    
    compute_final = True
    EPSI_CLUSTER = 25
    EPSI_ALPHASHAPE = 25
    USE_Z_ERASE = True
    IOU_THRESHOLD = 0.45
    print([compute_final, EPSI_CLUSTER, EPSI_ALPHASHAPE, USE_Z_ERASE, IOU_THRESHOLD])
    LAST_KK = 16
    onlyfiles1 = ["01_NI_Chil3-Cy3_Mki67-Cy5_02",
                 "04_IR5M_Chil3-Cy3_Serpine1-Cy5_06", 
                 "08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_010",
                 "11_NI_Cap-Cy3_Mki67-Cy5_002",
                 "07_CtrlNI_Pdgfra-Cy3_Serpine1-Cy5_002", 
                 "14_IR5M_Cap-Cy3_Serpine1-Cy5_011",
                 "03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_06",
                  "01_NI_Lamp3-Cy5_Pdgfra-Cy3_01",
                 "03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10"]
    
    path_to_json1 = [
        "/home/tom/Bureau/annotation/rna_association/01_NI_Chil3-Cy3_Mki67-Cy5_02_mip/target_files_v1/annotation.json", 
    "/home/tom/Bureau/annotation/rna_association/04_IR5M_Chil3-Cy3_Serpine1-Cy5_06_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_010_mip/target_files_v2/annotation.json",
     "/home/tom/Bureau/annotation/rna_association/11_NI_Cap-Cy3_Mki67-Cy5_002_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/07_CtrlNI_Pdgfra-Cy3_Serpine1-Cy5_002.json",
    "/home/tom/Bureau/annotation/rna_association/14_IR5M_Cap-Cy3_Serpine1-Cy5_011.json",
     "/home/tom/Bureau/annotation/rna_association/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06.json", 
    "/home/tom/Bureau/annotation/rna_association/01_NI_Lamp3-Cy5_Pdgfra-Cy3_01_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10_mip/target_files_v3/annotation.json"]
    
    
    onlyfiles2 = ["02_IR5M_Chil3-Cy3_Mki67-Cy5_05",
                 "03_NI_Chil3-Cy3_Serpine1-Cy5_01", 
                 "03_NI_Chil3-Cy3_Serpine1-Cy5_004",
                 "04_IR5M_Chil3-Cy3_Serpine1-Cy5_02",
                "12_IR5M_Cap-Cy3_Mki67-Cy5_005",
                "12_IR5M_Cap-Cy3_Mki67-Cy5_008",
                "12_IR5M_Cap-Cy3_Mki67-Cy5_009"
        ]
    
    path_to_json2 = [ "/home/tom/Bureau/annotation/rna_association/02_IR5M_Chil3-Cy3_Mki67-Cy5_05_mip/target_files_v1/annotation.json", 
    "/home/tom/Bureau/annotation/rna_association/03_NI_Chil3-Cy3_Serpine1-Cy5_01/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/03_NI_Chil3-Cy3_Serpine1-Cy5_004_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/04_IR5M_Chil3-Cy3_Serpine1-Cy5_02/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_005_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_008/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_009_mip/target_files_v2/annotation_c.json"
    ]
    
    
    path_to_json3 = path_to_json1 + path_to_json2 
    onlyfiles3  = onlyfiles1 + onlyfiles2
        

    list_res = []
    list_res_cluster = []
    list_res_cluster2 = []
    list_res_final = []    
    for kk in range(0,LAST_KK):
        dico_pred = {}
        dico_pred_cluster = {}
        dico_pred_cluster2 = {}
        dico_pred_final = {}
        print(kk)
        for index in range(len(onlyfiles3)):
            print(onlyfiles3[index])
            json_file = path_to_json3[index]
            print(json_file)
            with open(str(json_file), encoding='utf-8-sig') as fh:
                data_json = json.load(fh)
                label_Cy3_noise = geojson_to_label(data_json, img_size=(1040,1388),label = "Cy3_noise", binary_labeling=False)
                label_Cy5_noise = geojson_to_label(data_json, img_size=(1040,1388),label = "Cy5_noise",
                                                  binary_labeling=False)
        
                label_Cy3 = geojson_to_label(data_json, img_size=(1040,1388), label = "Cy3", binary_labeling=False)
                label_Cy5 = geojson_to_label(data_json, img_size=(1040,1388), label = "Cy5",
                                                  binary_labeling=False)
        
                binnary_label_Cy5_noise = label_Cy5_noise >= 1
                binnary_label_Cy3_noise = label_Cy3_noise >= 1
                binnary_label_Cy5 = label_Cy5 >= 1
                binnary_label_Cy3 = label_Cy3 >= 1
                
            spots_647 = np.load("/home/tom/Bureau/annotation/tiff_data2804/spots_af647/" + onlyfiles3[index] + ".npy")
            spots_568 = np.load("/home/tom/Bureau/annotation/tiff_data2804/spots_af568/" + onlyfiles3[index] + ".npy")
            rna_568 = tifffile.imread(path_rna_568+ "AF568_"+ onlyfiles3[index] + '.tiff')
            rna_647 = tifffile.imread(path_rna_647+ "AF647_"+ onlyfiles3[index] + '.tiff')

            new_spots_568, removed_spots_568, new_spots_647, removed_spots_647 = erase_overlapping_spot(spots_568, 
                                                                                spots_647, kk_568 = kk, kk_647 = kk,  use_z = USE_Z_ERASE)

            
            new_spots_568_cluster = erase_point_in_cluster_2Dalphashape(new_spots_568, removed_spots_568, eps=EPSI_ALPHASHAPE, 
                                                   min_samples = 4, min_cluster_size=4, xi=0.05,)
            new_spots_647_cluster = erase_point_in_cluster_2Dalphashape(new_spots_647, removed_spots_647, eps=EPSI_ALPHASHAPE, 
                                                 min_samples = 4, min_cluster_size=4, xi=0.05)
            
            if compute_final:
                
                labels_568 = computer_optics_cluster(new_spots_568_cluster, eps=EPSI_CLUSTER, min_samples = 4, min_cluster_size=4, xi=0.05)
                labels_647 = computer_optics_cluster(new_spots_647_cluster, eps=EPSI_CLUSTER, min_samples = 4, min_cluster_size=4, xi=0.05)
                
                img_dapi_mask = tifffile.imread("/home/tom/Bureau/annotation/tiff_data2804/predicted_mask_dapi/" +"dapi_maskdapi_"+ onlyfiles3[index] + '.tiff')

                
                nuclei_568_1, positive_cluster_568 = cluster_over_nuclei_3D_convex_hull(labels_568, 
                                                                                           np.array(new_spots_568_cluster),
                                                                                           img_dapi_mask, 
                                                                                            iou_threshold = IOU_THRESHOLD)
                nuclei_647_1, positive_cluster_647 =  cluster_over_nuclei_3D_convex_hull(labels_647,
                                                                                         np.array(new_spots_647_cluster), 
                                                                                             img_dapi_mask, iou_threshold = IOU_THRESHOLD)
                
                spots_568_after_clustering = []
                for s_index in range(len(new_spots_568_cluster)):
                    if labels_568[s_index] != -1:
                        spots_568_after_clustering.append(new_spots_568_cluster [s_index])
                        

                spots_647_after_clustering = []
                for s_index in range(len(new_spots_647_cluster)):
                    if labels_647[s_index] != -1:
                        spots_647_after_clustering.append(new_spots_647_cluster[s_index])
                
                list_cluster_568 = [p[0] for p in positive_cluster_568]
                list_cluster_647 = [p[0] for p in positive_cluster_647]
                
                final_point_set_568 = []
                for s_index in range(len(new_spots_568_cluster)):
                    if labels_568[s_index] in list_cluster_568:
                        final_point_set_568.append(new_spots_568_cluster[s_index])
                
                final_point_set_647 = []
                for s_index in range(len(new_spots_647_cluster)):
                    if labels_647[s_index] in list_cluster_647:
                        final_point_set_647.append(new_spots_647_cluster[s_index])

            def count_in_art(spots, labels):
                nb_of_spot_in_art = 0
                for s in spots:
                    if labels[s[1], s[2]]:
                       nb_of_spot_in_art += 1
                return nb_of_spot_in_art
            
            def count_in_real(spots, labels_real, label_art):
                nb_of_spot_in_label = 0
                for s in spots:
                    if labels_real[s[1], s[2]] and not label_art[s[1], s[2]]:
                        nb_of_spot_in_label += 1
                return nb_of_spot_in_label
            
            def get_in_art(spots, labels):
                spot_in_art = []
                for s in spots:
                    if labels[s[1], s[2]]:
                        spot_in_art.append(s)
                return spot_in_art
            
            def get_in_real(spots, labels_real, label_art):
                spot_in_label = []
                for s in spots:
                    if labels_real[s[1], s[2]] and not label_art[s[1], s[2]]:
                        spot_in_label.append(s)
                return np.array(spot_in_label)
            
            dico_pred[onlyfiles3[index]] = [np.array([count_in_art(spots_568, binnary_label_Cy3_noise),
                                             count_in_art(spots_647, binnary_label_Cy5_noise), 
                                             count_in_real(spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                             count_in_real(spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)]),
                                            
                                            np.array([count_in_art(new_spots_568, binnary_label_Cy3_noise),
                                            count_in_art(new_spots_647, binnary_label_Cy5_noise), 
                                             count_in_real(new_spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                             count_in_real(new_spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)])]
            
            dico_pred_cluster[onlyfiles3[index]] = [np.array([count_in_art(spots_568, binnary_label_Cy3_noise),
                                             count_in_art(spots_647, binnary_label_Cy5_noise), 
                                             count_in_real(spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                             count_in_real(spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)]),
                                            
                                            np.array([count_in_art(new_spots_568_cluster, binnary_label_Cy3_noise),
                                             count_in_art(new_spots_647_cluster, binnary_label_Cy5_noise), 
                                             count_in_real(new_spots_568_cluster, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                             count_in_real(new_spots_647_cluster, binnary_label_Cy5, binnary_label_Cy5_noise)])]
            if compute_final:
                
                
                dico_pred_cluster2[onlyfiles3[index]] = [np.array([count_in_art(spots_568, binnary_label_Cy3_noise),
                                                 count_in_art(spots_647, binnary_label_Cy5_noise), 
                                                 count_in_real(spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                 count_in_real(spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)]),
                                                
                                                np.array([count_in_art(spots_568_after_clustering, binnary_label_Cy3_noise),
                                                 count_in_art(spots_647_after_clustering, binnary_label_Cy5_noise), 
                                                 count_in_real(spots_568_after_clustering, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                count_in_real(spots_647_after_clustering, binnary_label_Cy5, binnary_label_Cy5_noise)])]
                
                
                dico_pred_final[onlyfiles3[index]] = [np.array([count_in_art(spots_568, binnary_label_Cy3_noise),
                                                 count_in_art(spots_647, binnary_label_Cy5_noise), 
                                                 count_in_real(spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                 count_in_real(spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)]),
                                                
                                                np.array([count_in_art(final_point_set_568, binnary_label_Cy3_noise),
                                                 count_in_art(final_point_set_647, binnary_label_Cy5_noise), 
                                                 count_in_real(final_point_set_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                count_in_real(final_point_set_647, binnary_label_Cy5, binnary_label_Cy5_noise)])]
        list_res.append(compute_pre_rec_mean(dico_pred))
        list_res_cluster.append(compute_pre_rec_mean(dico_pred_cluster))
        if compute_final:
            list_res_final.append(compute_pre_rec_mean(dico_pred_final, mean_bool = False))
            list_res_cluster2.append(compute_pre_rec_mean(dico_pred_cluster2, mean_bool = False))
    print([compute_final, EPSI_CLUSTER, EPSI_ALPHASHAPE, USE_Z_ERASE, IOU_THRESHOLD])
#%%   
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    fig.suptitle( "Cap-Cy3_Mki67-Cy5_", fontsize = 5)
    last_ind =18
    ax.plot(list(range(0, 25)), [np.mean(l[0]) for l in list_res], label = "cy3 recall ")

    ax.plot(list(range(0, 25)), [np.mean(l[1]) for l in list_res], label = "cy3 precision")
    ax.plot(list(range(0, 25)), [np.mean(l[4][4]/l[4][3]) for l  in list_res], label = "cy3 true negative rate ")
    ax.plot(list(range(0, 25)), [np.mean(l[0]) for l in list_res_cluster], label = "cy3 recall with α-shape suppression")
    ax.plot(list(range(0, 25)), [np.mean(l[1]) for l in list_res_cluster], label = "cy3 precision with α-shape suppression")
    ax.plot(list(range(0, 25)), [np.mean(l[0]) for l in list_res_cluster2], label = "cy3 recall after clustering")
    ax.plot(list(range(0, 25)), [np.mean(l[1]) for l in list_res_cluster2], label = "cy3 precision after clustering")
    ax.plot(list(range(0, 25)), [np.mean(l[0]) for l in list_res_final], label = "cy3 recall with final")
    ax.plot(list(range(0, 25)), [np.mean(l[1]) for l in list_res_final], label = "cy3 precision with final")
    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    plt.show() 
    
    
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    fig.suptitle( "Cap-Cy3_Mki67-Cy5", fontsize = 5)

    ax.plot(list(range(0, 25)), [np.mean(l[2]) for l in list_res], label = "cy5 recall")

    ax.plot(list(range(0, 25)), [np.mean(l[3]) for l in list_res], label = "cy5 precision")
    ax.plot(list(range(0, 25)), [np.mean(l[4][4]/l[4][3]) for l  in list_res], label = "cy3 true negative rate ")
    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    plt.show() 
    ax.plot(list(range(0, 25)), [np.mean(l[2]) for l in list_res_cluster], label = "cy5 recall with α-shape suppression")
    ax.plot(list(range(0, 25)), [np.mean(l[3]) for l in list_res_cluster], label = "cy5 precision with α-shape suppression")
    
    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    plt.show() 
    ax.plot(list(range(0, 25)), [np.mean(l[2]) for l in list_res_cluster2], label = "cy5 recall after clustering")
    ax.plot(list(range(0, 25)), [np.mean(l[3]) for l in list_res_cluster2], label = "cy5 precision after clustering")
    ax.plot(list(range(0, 25)), [np.mean(l[2]) for l in list_res_final], label = "cy5 recall with final")
    ax.plot(list(range(0, 25)), [np.mean(l[3]) for l in list_res_final], label = "cy5 precision  final")
    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    
    plt.show() 
    
#%%
    index_chill_serpine = [False,  True,  False,  False,  False,  False,  False,  False,  
                           False,  False,  True,  True,  True,  False,  False, False]
    index_extr = np.array([True] * 16)[6:9]
    nbkk =16
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    fig.suptitle( "", fontsize = 5)
    ax.plot(list(range(0, nbkk)), [np.mean(l[0][6:9]) for l in list_res[:nbkk]],
            label = "cy3 recall", color = "blue")
    ax.plot(list(range(0, nbkk)), [np.mean(l[1][6:9]) for l in list_res[:nbkk]], 
            label = "cy3 precision", color = "orange") 
    ax.plot(list(range(0, nbkk)), [np.mean(l[0][6:9]) for l in list_res_cluster],
            label = "cy3 recall after alphashape suppression",color ='red')
    ax.plot(list(range(0, nbkk)), [np.mean(l[1][6:9]) for l in list_res_cluster],
            label = "cy3 precision after alphashape suppression",  color ='purple')
    #ax.plot(list(range(0, nbkk)), [np.mean(l[0][6:9]) for l in list_res_final], label = "cy3 recall with final")
    #ax.plot(list(range(0, nbkk)), [np.mean(l[1][6:9]) for l in list_res_final], label = "cy3 precision with final") 
    #ax.plot(list(range(0, nbkk)), [np.mean(l[4][4][6:9]/l[4][3][6:9]) for l in list_res_final],  label = "cy3 true negative with final")

    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    plt.show() 
    
    
    
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    fig.suptitle( "", fontsize = 5)
    ax.plot(list(range(0, nbkk)), [np.mean(l[2][6:9]) for l in list_res[:nbkk]], label = "cy5 recall " , color = "blue")
    ax.plot(list(range(0, nbkk)), [np.mean(l[3][6:9]) for l in list_res[:nbkk]], label = "cy5 precision",  color = "orange")
    ax.plot(list(range(0, nbkk)), [np.mean(l[2][6:9]) for l in list_res_cluster], label = "cy5 recall after alphashape suppression", 
            color ='red')
    ax.plot(list(range(0, nbkk)), [np.mean(l[3][6:9]) for l in list_res_cluster], 
            label = "cy5 precision after alphashape suppression", color ='purple')
    #ax.plot(list(range(0, nbkk)), [np.mean(l[2]) for l in list_res_final], label = "cy5 recall with final")
    #ax.plot(list(range(0, nbkk)), [np.mean(l[3]) for l in list_res_final], label = "cy5 precision with final") 
    #ax.plot(list(range(0, nbkk)), [np.mean(l[5][4][l[5][3] != 0]/l[5][3][l[5][3] != 0] ) for l in list_res_final], label = "cy5 true negative with final")

    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    plt.show() 

#%% True negative rate plot
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    fig.suptitle( "Cy5 Lamp3, Pdgra", fontsize = 5)
    ax.plot(list(range(0, 8)), [np.mean(l[4][4]/l[4][3]) for l  in list_res], label = "cy3 true negative rate ")
    ax.plot(list(range(0, 8)), [np.mean(l[4][4]/l[4][3]) for l  in list_res_cluster], label = "cy3 true negativewith α-shape suppression")
    ax.plot(list(range(0, 8)), [np.mean(l[4][4]/l[4][3]) for l in list_res_final], label = "cy3 true negative with final")
    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    
    plt.show() 
    
    
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    fig.suptitle( "Cy5 Lamp3, Pdgra", fontsize = 5)

    ax.plot(list(range(0, 8)), [np.mean(l[5][4][:2]/l[5][3][:2]) for l in list_res], label = "cy5 tnr")
    ax.plot(list(range(0, 8)), [np.mean(l[5][4][:2]/l[5][3][:2]) for l  in list_res_cluster], label = "cy5 tnr with α-shape suppression")
    ax.plot(list(range(0, 8)), [np.mean(l[5][4][:2]/l[5][3][:2]) for  l in list_res_final], label = "cy5 tnr with final")
    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    
    plt.show() 
    

    
    
    
    #%%                   
    
        
   ### grid search
if __name__ == "__main__":
    path_to_json3 = ["/home/tom/Bureau/annotation/rna_association/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06.json", 
    "/home/tom/Bureau/annotation/rna_association/01_NI_Lamp3-Cy5_Pdgfra-Cy3_01_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10_mip/target_files_v3/annotation.json"]
    
    onlyfiles3  = ["03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_06", "01_NI_Lamp3-Cy5_Pdgfra-Cy3_01",  "03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10"]
    dico_grid = {}
    for minpts in [4]:
        for epsi in [10, 15, 20, 30, 35]:
            print(onlyfiles3)
            print(path_to_json3)
            list_res = []
            list_res_cluster = []
            for kk in range(1, 7):
                dico_pred = {}
                dico_pred_cluster = {}
                for index in range(len(onlyfiles3)):
                    print(onlyfiles3[index])
                    json_file = path_to_json3[index]
                    print(json_file)
                    with open(str(json_file), encoding='utf-8-sig') as fh:
                        data_json = json.load(fh)
                        label_Cy3_noise = geojson_to_label(data_json, img_size=(1040,1388),label = "Cy3_noise", binary_labeling=False)
                        label_Cy5_noise = geojson_to_label(data_json, img_size=(1040,1388),label = "Cy5_noise",
                                                          binary_labeling=False)
                
                        label_Cy3 = geojson_to_label(data_json, img_size=(1040,1388), label = "Cy3", binary_labeling=False)
                        label_Cy5 = geojson_to_label(data_json, img_size=(1040,1388), label = "Cy5",
                                                          binary_labeling=False)
                
                        binnary_label_Cy5_noise = label_Cy5_noise >= 1
                        binnary_label_Cy3_noise = label_Cy3_noise >= 1
                        binnary_label_Cy5 = label_Cy5 >= 1
                        binnary_label_Cy3 = label_Cy3 >= 1
                        
                    spots_647 = np.load("/home/tom/Bureau/annotation/tiff_data2804/spots_af647/" + onlyfiles3[index] + ".npy")
                    spots_568 = np.load("/home/tom/Bureau/annotation/tiff_data2804/spots_af568/" + onlyfiles3[index] + ".npy")
        
                    new_spots_568, removed_spots_568, new_spots_647, removed_spots_647 = erase_overlapping_spot(spots_568, 
                                                                                                                spots_647, kk_568 = kk, kk_647 = kk)
                    new_spots_568_cluster = erase_point_in_cluster_2Dalphashape(new_spots_568, removed_spots_568, eps=epsi, 
                                                           min_samples = 5, min_cluster_size=minpts, xi=0.05)
                    new_spots_647_cluster = erase_point_in_cluster_2Dalphashape(new_spots_647, removed_spots_647, eps=epsi, 
                                                           min_samples = 5, min_cluster_size=minpts, xi=0.05)
                    
            



                    dico_pred[onlyfiles3[index]] = [np.array([count_in_art(spots_568, binnary_label_Cy3_noise),
                                                     count_in_art(spots_647, binnary_label_Cy5_noise), 
                                                     count_in_real(spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                     count_in_real(spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)]),
                                                    
                                                    np.array([count_in_art(new_spots_568, binnary_label_Cy3_noise),
                                                     count_in_art(new_spots_647, binnary_label_Cy5_noise), 
                                                     count_in_real(new_spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                     count_in_real(new_spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)])]
                    
                    dico_pred_cluster[onlyfiles3[index]] = [np.array([count_in_art(spots_568, binnary_label_Cy3_noise),
                                                     count_in_art(spots_647, binnary_label_Cy5_noise), 
                                                     count_in_real(spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                     count_in_real(spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)]),
                                                    
                                                    np.array([count_in_art(new_spots_568_cluster, binnary_label_Cy3_noise),
                                                     count_in_art(new_spots_647_cluster, binnary_label_Cy5_noise), 
                                                     count_in_real(new_spots_568_cluster, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                     count_in_real(new_spots_647_cluster, binnary_label_Cy5, binnary_label_Cy5_noise)])]
                list_res.append(compute_pre_rec_mean(dico_pred))
                list_res_cluster.append(compute_pre_rec_mean(dico_pred_cluster))
            
                         
            dico_grid[(minpts,epsi)] = list_res_cluster
#%%
if __name__ == "__main__":
    nb_end = 16
    for minpts in [4]:
       
    
       for epsi in [10, 15, 20, 25, 30, 35, 40]:
           fig , ax = plt.subplots(1,1, figsize = (10, 10))
           fig.suptitle( str(minpts) + " "+ str(epsi))
           ax.plot(list(range(nb_end)), [l[0] for l in list_res], label = "cy3 recall")
           ax.plot(list(range(nb_end)), [l[1] for l in list_res], label = "cy3 precision")
           ax.plot(list(range(nb_end)), [l[0] for l in dico_grid[(minpts,epsi)]], label = "cy3 recall cluster "+str(minpts) + " "+ str(epsi))
           ax.plot(list(range(nb_end)), [l[1] for l in dico_grid[(minpts,epsi)]], label = "cy3 precision cluster "+str(minpts) + " "+ str(epsi))
           ax.set_xlabel('nearest neighbor euclidian dist in 100nm')
           ax.legend()           
           plt.show() 
           fig , ax = plt.subplots(1,1, figsize = (10, 10))
           fig.suptitle( str(minpts) + " "+ str(epsi))

           ax.plot(list(range(nb_end)), [l[2] for l in list_res], label = "cy5 recall")
           ax.plot(list(range(nb_end)), [l[3] for l in list_res], label = "cy3 precision")
           ax.plot(list(range(nb_end)), [l[2] for l in dico_grid[(minpts,epsi)]], label = "cy5 recall cluster"+str(minpts) + " "+ str(epsi))
           ax.plot(list(range(nb_end)), [l[3] for l in dico_grid[(minpts,epsi)]], label = "cy5 precision cluster"+str(minpts) + " "+ str(epsi))
           ax.set_xlabel('nearest neighbor euclidian dist in 100nm')
           ax.legend()
           plt.show() 

    
    
    
    

#%%grid search sur les parametres du clustering.







if __name__ == "__main__":
    
    compute_final = True
    
        
    
    onlyfiles1 = ["01_NI_Chil3-Cy3_Mki67-Cy5_02",
                 "04_IR5M_Chil3-Cy3_Serpine1-Cy5_06", 
                 "08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_010",
                 "11_NI_Cap-Cy3_Mki67-Cy5_002",
                 "07_CtrlNI_Pdgfra-Cy3_Serpine1-Cy5_002", 
                 "14_IR5M_Cap-Cy3_Serpine1-Cy5_011",
                 "03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_06",
                  "01_NI_Lamp3-Cy5_Pdgfra-Cy3_01",
                 "03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10"]
    
    path_to_json1 = [
        "/home/tom/Bureau/annotation/rna_association/01_NI_Chil3-Cy3_Mki67-Cy5_02_mip/target_files_v1/annotation.json", 
    "/home/tom/Bureau/annotation/rna_association/04_IR5M_Chil3-Cy3_Serpine1-Cy5_06_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_010_mip/target_files_v2/annotation.json",
     "/home/tom/Bureau/annotation/rna_association/11_NI_Cap-Cy3_Mki67-Cy5_002_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/07_CtrlNI_Pdgfra-Cy3_Serpine1-Cy5_002.json",
    "/home/tom/Bureau/annotation/rna_association/14_IR5M_Cap-Cy3_Serpine1-Cy5_011.json",
     "/home/tom/Bureau/annotation/rna_association/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06/03_IRM_Lamp3_Cy3_Pdgfra-Cy5_06.json", 
    "/home/tom/Bureau/annotation/rna_association/01_NI_Lamp3-Cy5_Pdgfra-Cy3_01_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/03_IR5M_Lamp3-Cy5_Pdgfra-Cy3_10_mip/target_files_v3/annotation.json"]
    
    
    onlyfiles2 = ["02_IR5M_Chil3-Cy3_Mki67-Cy5_05",
                 "03_NI_Chil3-Cy3_Serpine1-Cy5_01", 
                 "03_NI_Chil3-Cy3_Serpine1-Cy5_004",
                 "04_IR5M_Chil3-Cy3_Serpine1-Cy5_02",
                "12_IR5M_Cap-Cy3_Mki67-Cy5_005",
                "12_IR5M_Cap-Cy3_Mki67-Cy5_008",
                "12_IR5M_Cap-Cy3_Mki67-Cy5_009"
        ]
    
    path_to_json2 = [
    "/home/tom/Bureau/annotation/rna_association/02_IR5M_Chil3-Cy3_Mki67-Cy5_05_mip/target_files_v1/annotation.json", 
    "/home/tom/Bureau/annotation/rna_association/03_NI_Chil3-Cy3_Serpine1-Cy5_01/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/03_NI_Chil3-Cy3_Serpine1-Cy5_004_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/04_IR5M_Chil3-Cy3_Serpine1-Cy5_02/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_005_mip/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_008/target_files_v1/annotation.json",
    "/home/tom/Bureau/annotation/rna_association/12_IR5M_Cap-Cy3_Mki67-Cy5_009_mip/target_files_v2/annotation_c.json"
    ]
    
    
    path_to_json3 = path_to_json1 + path_to_json2 
    onlyfiles3  = onlyfiles1 + onlyfiles2

    
    list_res = []
    list_res_cluster2 = []
    list_res_final = []    
    for minpts in [4]:
        for epsi in [15,  25, 35,  40, 50]:
            dico_pred = {}
            dico_pred_cluster2 = {}
            dico_pred_final = {}
            for index in range(len(onlyfiles3)):
                print(onlyfiles3[index])
                json_file = path_to_json3[index]
                print(json_file)
                with open(str(json_file), encoding='utf-8-sig') as fh:
                    data_json = json.load(fh)
                    label_Cy3_noise = geojson_to_label(data_json, img_size=(1040,1388),label = "Cy3_noise", binary_labeling=False)
                    label_Cy5_noise = geojson_to_label(data_json, img_size=(1040,1388),label = "Cy5_noise",
                                                      binary_labeling=False)
            
                    label_Cy3 = geojson_to_label(data_json, img_size=(1040,1388), label = "Cy3", binary_labeling=False)
                    label_Cy5 = geojson_to_label(data_json, img_size=(1040,1388), label = "Cy5",
                                                      binary_labeling=False)
            
                    binnary_label_Cy5_noise = label_Cy5_noise >= 1
                    binnary_label_Cy3_noise = label_Cy3_noise >= 1
                    binnary_label_Cy5 = label_Cy5 >= 1
                    binnary_label_Cy3 = label_Cy3 >= 1
                
                spots_568 = np.load("/home/tom/Bureau/annotation/tiff_data2804/spots_af568/" + onlyfiles3[index] + ".npy")
                spots_647 = np.load("/home/tom/Bureau/annotation/tiff_data2804/spots_af647/" + onlyfiles3[index] + ".npy")
    
                labels_568 = computer_optics_cluster(spots_568, eps=epsi, min_samples = minpts, 
                                                     min_cluster_size=minpts, xi=0.05)
                labels_647 = computer_optics_cluster(spots_647, eps=epsi, min_samples = minpts, 
                                                     min_cluster_size=minpts, xi=0.05)
                
                img_dapi_mask = tifffile.imread("/home/tom/Bureau/annotation/tiff_data2804/predicted_mask_dapi/" +
                                                "dapi_maskdapi_"+ onlyfiles3[index] + '.tiff')
    
                
                nuclei_568_1, positive_cluster_568 = cluster_over_nuclei_3D_convex_hull(labels_568, 
                                                                                           spots_568, img_dapi_mask, 
                                                                                            iou_threshold = 0.45)
                nuclei_647_1, positive_cluster_647 =  cluster_over_nuclei_3D_convex_hull(labels_647, spots_647, 
                                                                                             img_dapi_mask, iou_threshold = 0.45)
                
                spots_568_after_clustering = []
                for s_index in range(len(spots_568)):
                    if labels_568[s_index] != -1:
                        spots_568_after_clustering.append(spots_568[s_index])
                        

                spots_647_after_clustering = []
                for s_index in range(len(spots_647)):
                    if labels_647[s_index] != -1:
                        spots_647_after_clustering.append(spots_647[s_index])

               
                list_cluster_568 = np.unique([p[0] for p in positive_cluster_568])
                list_cluster_647 = np.unique([p[0] for p in positive_cluster_647])
                
                final_point_set_568 = []
                for s_index in range(len(spots_568)):
                    if labels_568[s_index] in list_cluster_568:
                        final_point_set_568.append(spots_568[s_index])
                
                final_point_set_647 = []
                for s_index in range(len(spots_647)):
                    if labels_647[s_index] in list_cluster_647:
                        final_point_set_647.append(spots_647[s_index])


                def count_in_art(spots, labels):
                    nb_of_spot_in_art = 0
                    for s in spots:
                        if labels[s[1], s[2]]:
                           nb_of_spot_in_art += 1
                    return nb_of_spot_in_art
                
                def count_in_real(spots, labels_real, label_art):
                    nb_of_spot_in_label = 0
                    for s in spots:
                        if labels_real[s[1], s[2]] and not label_art[s[1], s[2]]:
                            nb_of_spot_in_label += 1
                    return nb_of_spot_in_label
                    

                dico_pred_final[onlyfiles3[index]] = [np.array([count_in_art(spots_568, binnary_label_Cy3_noise),
                                                 count_in_art(spots_647, binnary_label_Cy5_noise), 
                                                 count_in_real(spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                 count_in_real(spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)]),
                                                
                                                np.array([count_in_art(final_point_set_568, binnary_label_Cy3_noise),
                                                 count_in_art(final_point_set_647, binnary_label_Cy5_noise), 
                                                 count_in_real(final_point_set_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                 count_in_real(final_point_set_647, binnary_label_Cy5, binnary_label_Cy5_noise)])]
                
                dico_pred_cluster2[onlyfiles3[index]] = [np.array([count_in_art(spots_568, binnary_label_Cy3_noise),
                                                 count_in_art(spots_647, binnary_label_Cy5_noise), 
                                                 count_in_real(spots_568, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                 count_in_real(spots_647, binnary_label_Cy5, binnary_label_Cy5_noise)]),
                                                
                                                np.array([count_in_art(spots_568_after_clustering, binnary_label_Cy3_noise),
                                                 count_in_art(spots_647_after_clustering, binnary_label_Cy5_noise), 
                                                 count_in_real(spots_568_after_clustering, binnary_label_Cy3, binnary_label_Cy3_noise), 
                                                count_in_real(spots_647_after_clustering, binnary_label_Cy5, binnary_label_Cy5_noise)])]


            list_res_final.append(compute_pre_rec_mean(dico_pred_final,  mean_bool = False))
            list_res_cluster2.append(compute_pre_rec_mean(dico_pred_cluster2,   mean_bool = False))


    
            
            
            
    #%%     
            
            
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    fig.suptitle( "Cy3 Lamp3, Pdgra", fontsize = 5)
    
    ax.plot( [15, 20, 25, 30, 35,  40], [np.mean(l[0]) for l in list_res_final],
            label = "cy3 recall with final")
    ax.plot( [15, 20, 25, 30, 35,  40], [np.mean(l[1]) for l in list_res_final], 
            label = "cy3 precision with final")
    ax.plot( [15, 20, 25, 30, 35,  40], [np.mean(l[4][4]/l[4][3]) for l in list_res_final], label = "cy3 true negative with final")
    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    
    plt.show() 
    
    
    
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    fig.suptitle( "Cy3 Lamp3, Pdgra", fontsize = 5)
    
    ax.plot( [15, 20, 25, 30, 35,  40], [np.mean(l[2]) for l in list_res_final], label = "cy5 recall with final")
    ax.plot( [15, 20, 25, 30, 35,  40], [np.mean(l[3]) for l in list_res_final], label = "cy5 precision with final")
    ax.plot( [15, 20, 25, 30, 35,  40], [np.mean(l[5][4][:2]/l[5][3][:2]) for  l in list_res_final], label = "cy5 tnr with final")
    
    ax.set_xlabel('nearest neighbor euclidian dist (in 100nm)')
    ax.legend()
    
    plt.show() 
    
    
    
    #%%     
            
            
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    #fig.suptitle( "Cy3 Lamp3, Pdgra", fontsize = 5)
    
    ax.plot([15,  25, 35,  40, 50], [np.mean(l[0]) for l in list_res_final], label = "cy3 recall with final", color = "blue")
    ax.plot([15,  25, 35,  40, 50], [np.mean(l[1]) for l in list_res_final], label = "cy3 precision with final", color = "orange")
    ax.plot([15,  25, 35,  40, 50], [np.mean(l[4][4]/l[4][3]) for l in list_res_final], label = "cy3 true positive rate", color = "green")
    
    ax.set_xlabel('Espilon parameter of the DBSCAN', fontsize = 15)
    ax.legend()
    
    plt.show() 
    
    
    
    fig , ax = plt.subplots(1,1, figsize = (10, 10))
    
    ax.plot([15,  25, 35,  40, 50], [np.mean(l[2]) for l in list_res_final], label = "cy5 recall with final", color = "blue")
    ax.plot([15,  25, 35,  40, 50], [np.mean(l[3]) for l in list_res_final], label = "cy5 precision with final", color = "orange")
    ax.plot([15,  25, 35,  40, 50], [np.mean(l[5][4][l[5][3] != 0]/l[5][3][l[5][3] != 0]) for l in list_res_final], 
            label = "cy3 true positive rate", color = "green")

    ax.set_xlabel('Espilon parameter of the DBSCAN', fontsize = 15)
    ax.legend()
    
    plt.show() 




