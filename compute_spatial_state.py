#!/usr/bin/env python3

# -*- coding: utf-8 -*-
# %%
########################################
# File with fonction to generate excel describing the spatial state of cell type
##########################################

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage  # %%
from scipy import ndimage as ndi
from skimage.segmentation import watershed

from spots.post_processing import erase_solitary
from utils.utils import get_dye


def compute_average_nuclei_size(img_dapi_mask, positive_nuclei,
                                scale_z=300, scale_xy=103,
                                in_micrometre=True):
    """

    Args:
        img_dapi_mask ():
        positive_nuclei (list of nuclei positive to the probe):
        scale_z (float): in nanometer
        scale_xy (float): in nanometer
        in_micrometre (bool): if true convert the volume in micrometre

    Returns:
        return the mean estimated volume of the "positive" nuclei
    """

    total_volume = 0
    for nuc in positive_nuclei:
        total_volume += np.sum((img_dapi_mask == nuc).astype(np.int64)) * np.int64(scale_z) * np.int64(
            scale_xy) * np.int64(scale_xy)
    if in_micrometre:
        total_volume = total_volume * (10 ** -9)
    return total_volume / len(positive_nuclei) if len(positive_nuclei) > 0 else "Not define"


def get_experiment_name(sample_name):
    if "NI" in sample_name:
        return sample_name[:5]
    else:
        return sample_name[:7]


def count_positive_cell(dico_stat, image_name,
                        dye, exclude_impossible_solution=True):
    """

    Args:
        dico_stat (dict): dictionary computed in main_cluster.py
        image_name (str): probe name
        dye (str):
        exclude_impossible_solution (bool): parameter to not take into account nuclei that are positive to two
    incompatible cell type

    Returns:

    """

    cell_cy3 = [c[3] for c in dico_stat[image_name][5]]
    cell_cy5 = [c[3] for c in dico_stat[image_name][6]]
    print(type(image_name))
    print(image_name)
    if exclude_impossible_solution:
        if "Serpine1" in image_name or 'Mki67' in image_name: # cell state probe
            if dye == "Cy3":
                return len(cell_cy3), cell_cy3
            elif dye == "Cy5":
                return len(cell_cy5), cell_cy5
            else:
                raise Exception("Dye not detected")
        elif all(w in image_name for w in ['Pecam1', "Apln"]) or all(
                w in image_name for w in ['Pecam1', "Ptprb"]) or all(w in image_name for w in ['Hhip', "Pdgfra"]):
            if dye == "Cy3":
                return len(cell_cy3), cell_cy3
            elif dye == "Cy5":
                return len(cell_cy5), cell_cy5
            else:
                raise Exception("Dye not detected")
        else:
            if dye == "Cy3":
                return len([cell for cell in cell_cy3 if cell not in cell_cy5]), [cell for cell in cell_cy3 if
                                                                                  cell not in cell_cy5]
            elif dye == "Cy5":
                return len([cell for cell in cell_cy5 if cell not in cell_cy3]), [cell for cell in cell_cy5 if
                                                                                  cell not in cell_cy3]
            else:
                raise Exception("Dye not detected")
    else:
        if dye == "Cy3":
            return len(cell_cy3), cell_cy3
        elif dye == "Cy5":
            return len(cell_cy3), cell_cy5
        else:
            raise Exception("Dye not detected")


def compute_average_size(l_d, nuclei=None):
    """

    Args:
        l_d (list):  is the list of list [cluster number, overlapp, cluster volume, nuclei]
        nuclei (list): list of int of positive nuclei

    Returns:
        mean estimated point cloud size
    """

    if nuclei is not None and len(nuclei) == 0:
        return "not defined"
    if len(l_d) == 0:
        return "not defined"
    try:
        unique_clusters, frequency_cluster = np.unique(np.array(l_d)[:, 0], return_counts=True)
    except Exception as e:
        print(l_d)
        raise e
    dico_int = {}
    if nuclei is not None:
        l_d = [lll for lll in l_d if lll[3] in nuclei]
    for tup in l_d:
        dico_int[tup[0]] = []
    for tup in l_d:
        dico_int[tup[0]].append(tup[2] * (10 ** (-9)) / frequency_cluster[
            unique_clusters == tup[0]])  # add the point cloud volume of each cluster divide by the number of cell
    total_sum = 0
    for k in dico_int.keys():
        total_sum += np.sum(dico_int[k])
    if len(l_d) > 0:
        return total_sum / len(l_d)
    return "not defined"


def compute_average_size_precise(l_d, nuclei=None, percent_coverage=0.8):

    """
        take only into account well define point cloud where only one nucleus is present on the point cloud
    Args:
        l_d (list):  is the list of list (cluster number, overlapp, cluster volume, nuclei)
        nuclei (list): list of int of positive nuclei
        percent_coverage (float): minimal percentage coverage of the pint cloud over the nucleus
    Returns:
        mean estimated point cloud size
    """

    if len(l_d) == 0:
        return None, "not defined"
    unique_clusters, frequency_cluster = np.unique(np.array(l_d)[:, 0], return_counts=True)
    if nuclei is not None:
        l_d = [lll for lll in l_d if lll[3] in nuclei]
    cluster_freq_1 = unique_clusters[frequency_cluster == 1]
    if len(l_d) == 0:
        return None, "not defined"
    l_d = [lll for lll in l_d if lll[0] in cluster_freq_1 and lll[1] > percent_coverage]
    if len(l_d) == 0:
        return None, "not defined"
    return np.sum([vol[2] * (10 ** -9) for vol in l_d]) / len(l_d), len(l_d)




# %%

def generate_exels_one_cell(list_folder="",
                            gene_smfish=('Cap', 'aCap', 'CEC', 'acap'),
                            path_to_take="/media/tom/Transcend/image_lustra0605/Images_Hugo/",
                            path_save="/media/tom/Transcend/image_lustra0605/exelsonecell/",
                            dico_stat_name="dico_seg1005.npy",
                            compute_nuclei_size=False,
                            scale_z=300,
                            scale_xy=103
                            ):
    """
    Parameters
    ----------
    list_param
    Returns
    -------
    """

    columns = [
        "folder_name",
        "experiment",
        "image_name",
        "gene",
        "dye",
        "nb_nuclei",
        "nb_positive",
        "average_point_cloud_size",
        "average_nuclei_size",
        "average_point cloud " + str(gene_smfish[0]) + " with point clouds containing only one nucleus",
        'number of nuclei used for ' + str(gene_smfish[0]) + " estimation",
    ]
    dataframe_per_files_NI = pd.DataFrame(columns=columns)

    dataframe_per_files_IR1M = pd.DataFrame(columns=columns)

    dataframe_per_files_IR2M = pd.DataFrame(columns=columns)

    dataframe_per_files_IR3M = pd.DataFrame(columns=columns)

    dataframe_per_files_IR4M = pd.DataFrame(columns=columns)

    dataframe_per_files_IR5M = pd.DataFrame(columns=columns)

    dataframe_per_files_other_IRM = pd.DataFrame(columns=columns)

    for folder_name in list_folder:
        print(folder_name)
        path_output_segmentaton = path_to_take + folder_name + "tiff_data/" + "predicted_mask_dapi/"
        dico_stat = np.load(path_to_take + folder_name + dico_stat_name, allow_pickle=True).item()
        onlyfiles = [list(dico_stat.keys())]
        print(onlyfiles)
        sorted_name = np.sort(list(dico_stat.keys()))
        for image_name in sorted_name:
            if not any(word in image_name for word in gene_smfish):
                continue
            print(image_name)
            experience_name = get_experiment_name(image_name)
            dye = get_dye(gene_smfish, image_name)
            nb_positive, positive_nuclei = count_positive_cell(dico_stat, image_name, dye)
            average_point_cloud_size = compute_average_size(dico_stat[image_name][5],
                                                            positive_nuclei) if dye == 'Cy3' else compute_average_size(
                dico_stat[image_name][6],
                positive_nuclei)
            average_nuclei_size = None
            score_limit_list_over_neibors = []
            average_precise_point_cloud_size_type = None
            nb_c_type = None
            # shape index computation
            if compute_nuclei_size:
                img_dapi_mask = tifffile.imread(path_output_segmentaton + "dapi_maskdapi_" + image_name)
                img_dapi_mask = erase_solitary(img_dapi_mask)
                average_nuclei_size = compute_average_nuclei_size(img_dapi_mask, positive_nuclei, scale_z=scale_z,
                                                                  scale_xy=scale_xy, in_micrometre=True)

            dico_input_pd = {
                "folder_name": Path(folder_name).parts[-2:],
                "experiment": get_experiment_name(image_name),
                "image_name": image_name,
                "gene": gene_smfish,
                "dye": dye,
                "nb_nuclei": dico_stat[image_name][0],
                "nb_positive": nb_positive,
                "average_point_cloud_size": average_point_cloud_size,
                "average_nuclei_size": average_nuclei_size,
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
                      pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)],
                                   columns=dataframe_per_files_NI.columns),
                      dataframe_per_files_other_IRM,
                      pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)],
                                   columns=dataframe_per_files_NI.columns),
                      dataframe_per_files_IR5M]
            result = pd.concat(frames)

            if not os.path.exists(path_save):
                os.mkdir(path_save)
            print(result)
            result.to_pickle(path_save + str(gene_smfish[0]) + ".pkl")
            result.to_excel(path_save + str(gene_smfish[0]) + '.xls')
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
    result.to_pickle(path_save + str(gene_smfish[0]) + ".pkl")
    result.to_excel(path_save + str(gene_smfish[0]) + '.xls')
    print("save")
    return result


# %%
def generate_exels_two_probes(list_folder,
                                   gene_1,
                                   gene_2,
                                   path_to_take,
                                   path_save,
                                   dico_stat_name,
                                   scale_z=300,
                                   scale_xy=103,
                                   compute_nuclei_size=False):
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
        "average_nuclei_size  (μm3)" + str(gene_1[0]),

        "average_point_cloud_size (μm3)" + str(gene_2[0]),
        "average_nuclei_size (μm3)" + str(gene_2[0]),

        "average_point_cloud_size (μm3) positive_to_both",
        "average_nuclei_size (μm3) positive_to_both",

    ]
    dataframe_per_files_NI = pd.DataFrame(columns=columns)
    dataframe_per_files_IR5M = pd.DataFrame(columns=columns)
    dataframe_per_files_other_IRM = pd.DataFrame(columns=columns)
    for folder_name in list_folder:
        dico_stat = np.load(path_to_take + folder_name + dico_stat_name, allow_pickle=True).item()
        sorted_name = np.sort(list(dico_stat.keys()))
        for image_name in sorted_name:
            if not any(word in image_name for word in gene_1):
                continue
            if not any(word in image_name for word in gene_2):
                continue
            print(image_name)
            # gene_1
            dye_type = get_dye(gene_1, image_name)
            nb_positive_gene1, positive_nuclei_type = count_positive_cell(dico_stat, image_name, dye_type)
            average_point_cloud_size_type = compute_average_size(dico_stat[image_name][5],
                                                                 positive_nuclei_type) if dye_type == 'Cy3' else compute_average_size(
                dico_stat[image_name][6],
                positive_nuclei_type)

            if compute_nuclei_size:
                path_output_segmentaton = path_to_take + folder_name + "tiff_data/" + "predicted_mask_dapi/"
                img_dapi_mask = tifffile.imread(path_output_segmentaton + "dapi_maskdapi_" + image_name)
                img_dapi_mask = erase_solitary(img_dapi_mask)
                average_nuclei_size_gene_1 = compute_average_nuclei_size(img_dapi_mask, positive_nuclei_type,
                                                                       scale_z=scale_z, scale_xy=scale_xy)
            else:
                average_nuclei_size_gene_1 = None


            # gene_2
            dye_gene = get_dye(gene_2, image_name)
            nb_positive_gene2, positive_nuclei_state = count_positive_cell(dico_stat, image_name, dye_gene)

            average_point_cloud_size_state = compute_average_size(dico_stat[image_name][5],
                                                                  positive_nuclei_state) if dye_gene == 'Cy3' else compute_average_size(
                dico_stat[image_name][6], positive_nuclei_state)

            if compute_nuclei_size:
                average_nuclei_size_gene_2 = compute_average_nuclei_size(img_dapi_mask, positive_nuclei_state,
                                                                        scale_z=scale_z, scale_xy=scale_xy)
            else:
                average_nuclei_size_gene_2 = None

            # gene_both
            positive_both = list(set(positive_nuclei_state) & set(positive_nuclei_type))
            average_point_cloud_size_both = compute_average_size(dico_stat[image_name][5],
                                                                 positive_both)

            if compute_nuclei_size:
                average_nuclei_size_both = compute_average_nuclei_size(img_dapi_mask, positive_both, scale_z=scale_z,
                                                                       scale_xy=scale_xy)

            dico_input_pd = {

                "folder_name": Path(folder_name).parts[-2:],
                "experiment": get_experiment_name(image_name),
                "image_name": image_name,
                "probe1": gene_1[0],
                "probe2 ": gene_2[0],

                "nb_nuclei": dico_stat[image_name][0],
                "probe1 only": nb_positive_gene1,
                "probe2 only": nb_positive_gene2,
                "nb_positive_both": len(positive_both),

                "average_point_cloud_size (μm3) " + str(gene_1[0]): average_point_cloud_size_type,
                "average_nuclei_size  (μm3)" + str(gene_1[0]):  average_nuclei_size_gene_1,

                "average_point_cloud_size (μm3)" + str(gene_2[0]): average_point_cloud_size_state,
                "average_nuclei_size (μm3)" + str(gene_2[0]): average_nuclei_size_gene_2,

                "average_point_cloud_size (μm3) positive_to_both": average_point_cloud_size_both,
                "average_nuclei_size (μm3) positive_to_both": average_nuclei_size_both,

            }

            experience_name = get_experiment_name(image_name)
            if "NI" in experience_name:
                dataframe_per_files_NI.loc[len(dataframe_per_files_NI)] = dico_input_pd
            elif "IR5M" in experience_name:
                dataframe_per_files_IR5M.loc[len(dataframe_per_files_IR5M)] = dico_input_pd
            else:
                dataframe_per_files_other_IRM.loc[len(dataframe_per_files_other_IRM)] = dico_input_pd
            frames = [dataframe_per_files_NI,
                      pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)],
                                   columns=dataframe_per_files_NI.columns),
                      dataframe_per_files_other_IRM,
                      pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)],
                                   columns=dataframe_per_files_NI.columns),
                      dataframe_per_files_IR5M]
            result = pd.concat(frames)
            if not os.path.exists(path_save):
                os.mkdir(path_save)
            result.to_pickle(path_save + str(gene_1[0]) + "_" + gene_2[0] + ".pkl")
            result.to_excel(path_save + str(gene_1[0]) + "_" + gene_2[0] + '.xls')
    frames = [dataframe_per_files_NI,
              pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
              dataframe_per_files_other_IRM,
              pd.DataFrame([[''] * len(dataframe_per_files_NI.columns)], columns=dataframe_per_files_NI.columns),
              dataframe_per_files_IR5M]
    result = pd.concat(frames)
    if not os.path.exists(path_save):
        os.mkdir(path_save)
    result.to_pickle(path_save + str(gene_1[0]) + "_" + gene_2[0] + ".pkl")
    result.to_excel(path_save + str(gene_1[0]) + "_" + gene_2[0] + '.xls')
    print("save")
    return result
