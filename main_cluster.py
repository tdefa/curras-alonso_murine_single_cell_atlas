# -*- coding: utf-8 -*-

"""
Created on Wed Mar 31 11:27:43 2021
@author: thomas
"""
#%%

import argparse
import datetime
import gc
import json
import os
import re
import time
import warnings
from os import listdir
from os.path import isfile, join

import numpy as np
import tifffile
from cellpose import models
from matplotlib import pyplot as plt
from tqdm import tqdm

from run_seg import segment_nuclei
from spots.erase_overlapping_spot import (erase_overlapping_spot,
                                          )
from spots.plot import (
    mask_image_to_rgb,
    mask_image_to_rgb2D_from_list_green_cy3_red_cy5_both_blue_grey,
    mask_image_to_rgb2D_from_list_orange_cy3_other_grey,
    mask_image_to_rgb2D_from_list_orange_cy5_other_grey)
from spots.post_processing import erase_solitary
from spots.spot_detection import (cluster_over_nuclei_3D_convex_hull,
                                  computer_optics_cluster,
                                  spot_detection_for_clustering)
from utils.czi_to_tiff import preprare_tiff

warnings.filterwarnings("ignore")


def main(args):
    list_folder = args.list_folder
    dico_cy3 = json.loads(args.manual_threshold_cy3)
    dico_cy5 = json.loads(args.manual_threshold_cy5)
    print(f"manual threshold {dico_cy3} {dico_cy5} ")
    # default probe parameters
    dico_param_probes = {"Lamp3": (32, 0.42),
                         "Pdgfra": (35, 0.42),
                         "Chil3": (20, 0.55),
                         'Cap': (35, 0.30),
                         'aCap': (35, 0.30),
                         'acap': (35, 0.30),
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
    # modify default probe parameters or add new probes
    for new_probe in args.new_probe:
        dico_param_probes[new_probe[0]] = (int(new_probe[1]), float(new_probe[2]))
    print(f'list of probe names you can use {list(dico_param_probes.keys())}')

    # create the folder structure to save all intermediate and final results and plots
    print(f"input folder {list_folder}")
    for folder_index in range(len(list_folder)):
        print(list_folder[folder_index])
        folder = list_folder[folder_index]
        path_to_save_fig = args.path_to_czi_folder + "figure/"
        path_save_res_classif = args.path_to_czi_folder + "res_classif/"
        path_to_czi_folder_c = args.path_to_czi_folder + list_folder[folder_index]
        path_to_project_c = args.path_to_czi_folder + list_folder[folder_index]

        if not os.path.exists(path_to_czi_folder_c + "tiff_data/"):
            os.mkdir(path_to_czi_folder_c + "tiff_data/")
        if not os.path.exists(path_to_czi_folder_c + "tiff_data/" + "dapi/"):
            os.mkdir(path_to_czi_folder_c + "tiff_data/" + "dapi/")
        if not os.path.exists(path_to_czi_folder_c + "tiff_data/" + "af568/"):
            os.mkdir(path_to_czi_folder_c + "tiff_data/" + "af568/")
        if not os.path.exists(path_to_czi_folder_c + "tiff_data/" + "af647/"):
            os.mkdir(path_to_czi_folder_c + "tiff_data/" + "af647/")

        path_to_czi = path_to_czi_folder_c
        path_to_dapi = path_to_czi_folder_c + "tiff_data/" + "dapi/"
        path_to_af647 = path_to_czi_folder_c + "tiff_data/" + "af647/"
        path_to_af568 = path_to_czi_folder_c + "tiff_data/" + "af568/"
        path_output_segmentaton = path_to_czi_folder_c + "tiff_data/" + "predicted_mask_dapi/"

        #################
        # CZI TO tiff
        #################
        if args.prepare_czi:
            print('prepare czi')
            try:
                preprare_tiff(path_to_czi, path_to_dapi, path_to_af647, path_to_af568)
            except Exception as e:
                print(e)

        if not os.path.exists(path_output_segmentaton):
            os.mkdir(path_output_segmentaton)
        ###########
        # Nuclei segmentation
        ###########
        if args.segmentation:
            model = models.Cellpose(gpu=args.gpu, model_type='nuclei')
            # ##parameter
            dico_param = {}
            dico_param["diameter"] = args.diameter
            dico_param["gpu"] = bool(args.gpu)
            dico_param["flow_threshold"] = args.flow_threshold
            dico_param["do_3D"] = args.do_3D
            dico_param["mip"] = args.mip
            dico_param["projected_focused"] = False
            dico_param["stitch_threshold"] = args.stitch_threshold
            segment_nuclei(path_to_dapi, path_output_segmentaton, dico_param, model)

        ###########
        # check that all the nuclei are segmented
        ###########
        onlyfiles = [f for f in listdir(path_output_segmentaton) if
                     isfile(join(path_output_segmentaton, f)) and f[-1] == "f"]
        onlyfiles = [onlyfiles[i][14:] for i in range(len(onlyfiles))]
        assert len(onlyfiles) <= len([f for f in listdir(path_to_dapi) if isfile(join(path_to_dapi, f)) and f[-1] == "f"])



        ###########
        # Spot detection
        ###########


        if args.spot_detection:
            if not os.path.exists(path_to_project_c + "detected_spot_3d" + "/"):
                os.mkdir(path_to_project_c + "detected_spot_3d" + "/")
            dico_threshold = spot_detection_for_clustering(sigma=[1.25, 1.25, 1.25],
                                                           rna_path=[path_to_af568 + 'AF568_'],
                                                           path_output_segmentaton=path_output_segmentaton,
                                                           threshold_input=dico_cy3,
                                                           output_file=path_to_project_c + "detected_spot_3d" + "/", )
            np.save(path_to_project_c + 'dico_threshold_AF568.npy', dico_threshold)
            with open(path_to_project_c + 'dico_threshold_AF568.txt', 'w') as f:
                f.write(str(dico_threshold))

            print(dico_threshold)
            # different sigma than for af568 as the wave lenght is different
            dico_threshold = spot_detection_for_clustering(sigma=[1.35, 1.35, 1.35],
                                                           rna_path=[path_to_af647 + 'AF647_'],
                                                           path_output_segmentaton=path_output_segmentaton,
                                                           threshold_input=dico_cy5,
                                                           output_file=path_to_project_c + "detected_spot_3d" + "/", )
            np.save(path_to_project_c + 'dico_threshold_AF647.npy', dico_threshold)
            with open(path_to_project_c + 'dico_threshold_AF647.txt', 'w') as f:
                f.write(str(dico_threshold))

            print(dico_threshold)

        if args.classify:
            print("classify")
            try:  # load available result of available
                dico_stat = np.load(path_to_project_c + args.dict_name_save + '.npy', allow_pickle=True).item()
            except Exception as e:
                print(e)
                dico_stat = {}
            try:  # load available result of available
                dico_label_cluster = np.load(path_to_project_c + "dico_label_cluster" + args.dict_name_save + '.npy',
                                             allow_pickle=True).item()
            except Exception as e:
                print(e)
                dico_label_cluster = {}
            t = time.time()
            print(onlyfiles)
            for f in tqdm(onlyfiles[:]):
                print(f)
                print(list_folder[folder_index])
                print(time.time() - t)
                print(f[:-5])
                t = time.time()
                print(list(dico_label_cluster.keys()))

                args.epsi_cluster_cy3 = 'probe not reconize yet'
                args.epsi_cluster_cy5 = 'probe not reconize yet'

                # set clustering parameters
                for probe_name in dico_param_probes.keys():
                    if probe_name + '-Cy3' in f:
                        args.epsi_cluster_cy3 = dico_param_probes[probe_name][0]
                        args.overlapping_cy3 = dico_param_probes[probe_name][1]
                    elif probe_name + '-Cy5' in f:
                        args.epsi_cluster_cy5 = dico_param_probes[probe_name][0]
                        args.overlapping_cy5 = dico_param_probes[probe_name][1]
                    elif probe_name + 'Cy3' in f:
                        args.epsi_cluster_cy3 = dico_param_probes[probe_name][0]
                        args.overlapping_cy3 = dico_param_probes[probe_name][1]
                    elif probe_name + 'Cy5' in f:
                        args.epsi_cluster_cy5 = dico_param_probes[probe_name][0]
                        args.overlapping_cy5 = dico_param_probes[probe_name][1]
                    elif probe_name in f:
                        a = re.search(f'{probe_name}', f)
                        remaining_name = f[a.end():]
                        # if there are a probes after it means that probe_name is the first so it is Cy5
                        if any(word in remaining_name for word in list(dico_param_probes.keys())):
                            args.epsi_cluster_cy5 = dico_param_probes[probe_name][0]
                            args.overlapping_cy5 = dico_param_probes[probe_name][1]
                        # the second name is  Cy3
                        else:
                            args.epsi_cluster_cy3 = dico_param_probes[probe_name][0]
                            args.overlapping_cy3 = dico_param_probes[probe_name][1]
                if args.epsi_cluster_cy3 == 'probe not reconize yet':
                    raise Exception(f" probe not reconize for the file {f} in Cy3, set the argument new_probe ")
                if args.epsi_cluster_cy5 == 'probe not reconize yet':
                    raise Exception(f" probe not reconize for the file {f} in Cy5, set the argument new_probe ")
                assert args.epsi_cluster_cy3 != 'probe not reconize yet'
                assert args.epsi_cluster_cy5 != 'probe not reconize yet'
                ###########
                # load mask, remove too small nuclei segmentation artefact and load detected spots
                ###########
                print(f)
                img_dapi_mask = tifffile.imread(path_output_segmentaton + "dapi_maskdapi_" + f)
                img_dapi_mask = erase_solitary(img_dapi_mask)

                spots_568 = np.load(path_to_czi_folder_c + "detected_spot_3d" +
                                    "/" + "AF568_" + f[:-5] + 'array.npy')
                spots_647 = np.load(path_to_czi_folder_c + "detected_spot_3d" +
                                    "/" + "AF647_" + f[:-5] + 'array.npy')
                ######
                # erase possible artefact detection
                ######
                if args.remove_overlaping:  # erase overlapping detection
                    new_spots_568, removed_spots_568, new_spots_647, removed_spots_647 = erase_overlapping_spot(
                        spots_568,
                        spots_647,
                        kk_568=args.kn_568,
                        kk_647=args.kn_647,
                        scale=args.scale)
                    spots_568, spots_647 = np.array(new_spots_568), np.array(new_spots_647)

                ########
                # compute clustered dbscan
                #########
                if img_dapi_mask.ndim == 2:
                    spots_568 = np.array([[s[1], s[2]] for s in list(spots_568)])
                    spots_647 = np.array([[s[1], s[2]] for s in list(spots_647)])
                print(args.epsi_cluster_cy3)
                print(args.epsi_cluster_cy5)
                print(len(spots_568))
                print(len(spots_647))
                labels_568 = np.array([-1] * len(spots_568))
                labels_647 = np.array([-1] * len(spots_647))
                if type(args.epsi_cluster_cy3) == int:
                    labels_568 = computer_optics_cluster(spots_568, eps=args.epsi_cluster_cy3, min_samples=4,
                                                         min_cluster_size=4, xi=0.05)
                if type(args.epsi_cluster_cy5) == int:
                    labels_647 = computer_optics_cluster(spots_647, eps=args.epsi_cluster_cy5, min_samples=4,
                                                         min_cluster_size=4, xi=0.05)

                ##########
                # classify cell type
                #########
                if img_dapi_mask.ndim == 3:
                    nuclei_568_1, positive_cluster_568, negative_cluster_568 = cluster_over_nuclei_3D_convex_hull(
                        labels_568,
                        spots_568,
                        img_dapi_mask,
                        iou_threshold=args.overlapping_cy3)
                    nuclei_647_1, positive_cluster_647, negative_cluster_647 = cluster_over_nuclei_3D_convex_hull(
                        labels_647, spots_647,
                        img_dapi_mask, iou_threshold=args.overlapping_cy5)

                    nb_no_rna = len(np.unique(img_dapi_mask)) - len(set(nuclei_647_1).union(set(nuclei_568_1)))
                    nb_cy3 = len(set(nuclei_568_1) - set(nuclei_647_1))
                    nb_cy5 = len(set(nuclei_647_1) - set(nuclei_568_1))
                    nb_both = len(set(nuclei_647_1).intersection(set(nuclei_568_1)))
                    dico_stat[f] = [len(np.unique(img_dapi_mask)), nb_no_rna, nb_cy3, nb_cy5, nb_both,
                                    positive_cluster_568, positive_cluster_647, negative_cluster_568,
                                    negative_cluster_647]
                    # save spot but  without the overlapping spots
                    np.save(path_to_project_c + "dico_label_cluster" + args.dict_name_save, dico_label_cluster)
                    np.save(path_to_project_c + args.dict_name_save, dico_stat)

                dico_label_cluster[f] = [labels_568, labels_647, spots_568, spots_647]
                assert img_dapi_mask.ndim == 3

            np.save(path_to_project_c + "dico_label_cluster_final" + args.dict_name_save, dico_label_cluster)
            np.save(path_to_project_c + args.dict_name_save, dico_stat)

            if not os.path.exists(path_save_res_classif):
                os.mkdir(path_save_res_classif)
            np.save(path_save_res_classif + folder[:-2] + "dico_label_cluster_final" + args.dict_name_save,
                    dico_label_cluster)
            np.save(path_save_res_classif + folder[:-2]  + args.dict_name_save, dico_stat)
        #############
        # PLOTTING
        #############
        if args.save_plot:
            print("plotting")
            print(onlyfiles)
            print(onlyfiles)
            for f in onlyfiles:
                plt.close("all")
                print(f)
                if not os.path.exists(path_to_save_fig):
                    os.mkdir(path_to_save_fig)
                path_to_create_plot = ""
                for subfolder_plot in folder.split('/'):
                    path_to_create_plot += subfolder_plot +'/'
                    if not os.path.exists(path_to_save_fig + path_to_create_plot):
                        os.mkdir(path_to_save_fig + path_to_create_plot)
                # load analysis output
                dico_label_cluster = np.load(path_to_project_c + "dico_label_cluster_final" + args.dict_name_save + ".npy",
                                             allow_pickle=True).item()
                [labels_568, labels_647, spots_568, spots_647] = dico_label_cluster[f]
                dico_stat = np.load(path_to_project_c + args.dict_name_save + ".npy", allow_pickle=True).item()
                [nb_nuclei, nb_no_rna, nb_cy3, nb_cy5, nb_both, positive_cluster_568, positive_cluster_647, negative_cluster_568, negative_cluster_647] = dico_stat[f]
                nuclei_568_1 = [dico_stat[f][5][i][3] for i in range(len(dico_stat[f][5]))]
                nuclei_647_1 = [dico_stat[f][6][i][3] for i in range(len(dico_stat[f][6]))]

                img = tifffile.imread(path_to_dapi + "dapi_" + f)
                img_dapi_mask = tifffile.imread(path_output_segmentaton + "dapi_maskdapi_" + f)
                img_dapi_mask = erase_solitary(img_dapi_mask)

                #############
                # plot final classification
                #############

                if not os.path.exists(path_to_save_fig + folder + "classif/"):
                    os.mkdir(path_to_save_fig + folder + "classif/")
                fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                fig.suptitle(f + " grey (no rna) %s, Cy3+ green %s, Cy5+ red %s,  Both blue %s" % (str(nb_no_rna),
                                                                                                 str(nb_cy3),
                                                                                                 str(nb_cy5),
                                                                                                 str(nb_both)),
                             fontsize=20)

                m, green, yellow, blue, purple = mask_image_to_rgb2D_from_list_green_cy3_red_cy5_both_blue_grey(
                    np.amax(img, 0),
                    np.amax(img_dapi_mask, 0), nuclei_568_1, nuclei_647_1)
                ax.imshow(m)
                fig.savefig(path_to_save_fig + folder + "classif/green_cy3_red_cy5" + f[:-5])

                fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                fig.suptitle(f + "Cy3+ orange %s, other grey %s " % (str(nb_cy3),
                                                                    str(nb_cy5 + nb_no_rna + nb_both)), fontsize=20)

                m, green, yellow, blue, purple = mask_image_to_rgb2D_from_list_orange_cy3_other_grey(np.amax(img, 0),
                                                                                                     np.amax(
                                                                                                         img_dapi_mask,
                                                                                                         0),
                                                                                                     nuclei_568_1,
                                                                                                     nuclei_647_1)
                ax.imshow(m)
                fig.savefig(path_to_save_fig + folder + "classif/orange_cy3" + f[:-5])

                fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                fig.suptitle(f + "Cy5+ orange %s, other grey %s " % (str(nb_cy5),
                                                                    str(nb_cy3 + nb_no_rna + nb_both)), fontsize=20)

                m, green, yellow, blue, purple = mask_image_to_rgb2D_from_list_orange_cy5_other_grey(np.amax(img, 0),
                                                                                                     np.amax(
                                                                                                         img_dapi_mask,
                                                                                                         0),
                                                                                                     nuclei_568_1,
                                                                                                     nuclei_647_1)
                ax.imshow(m)
                fig.savefig(path_to_save_fig + folder + "classif/orange_cy5" + f[:-5])
                #####
                # plot segmentation
                #####
                if not os.path.exists(path_to_save_fig + folder + "segmentation/"):
                    os.mkdir(path_to_save_fig + folder + "segmentation/")

                fig, ax = plt.subplots(1, 1, figsize=(30, 20))

                m = mask_image_to_rgb(np.amax(img, 0), np.amax(img_dapi_mask, 0))
                ax.imshow(m)
                fig.savefig(path_to_save_fig + folder + "segmentation/" + f[:-5])
                #####
                # plot dapi_slpot ### Cy3 green cy5 orange
                #####
                if not os.path.exists(path_to_save_fig + folder + "dapi_spots/"):
                    os.mkdir(path_to_save_fig + folder + "dapi_spots/")

                # raw spot detection
                path_to_af568 = path_to_czi_folder_c + "tiff_data/" + "af568/"

                path_to_af647 = path_to_czi_folder_c + "tiff_data/" + "af647/"
                cy3_im = tifffile.imread(path_to_af568 + "AF568_" + f)
                cy5_im = tifffile.imread(path_to_af647 + "AF647_" + f)
                fig, ax = plt.subplots(2, 1, figsize=(30, 60))
                plt.title(f + " Cy3 fish green spots", fontsize=20)
                if img.ndim == 3:
                    ax[0].imshow(np.amax(cy3_im, 0))
                    ax[1].imshow(np.amax(cy3_im, 0))
                else:
                    ax[0].imshow(cy3_im)
                    ax[1].imshow(cy3_im)
                for s in spots_568:
                    ax[0].scatter(s[-1], s[-2], c='green', s=35)
                fig.savefig(path_to_save_fig + folder + "dapi_spots/smfish_cy3" + f[:-5])

                fig, ax = plt.subplots(2, 1, figsize=(30, 60))
                plt.title(f + " Cy5 fish green spots", fontsize=28)
                if img.ndim == 3:
                    ax[0].imshow(np.amax(cy5_im, 0))
                    ax[1].imshow(np.amax(cy5_im, 0))
                else:
                    ax[0].imshow(cy5_im)
                    ax[1].imshow(np.amax(cy5_im, 0))
                for s in spots_647:
                    ax[0].scatter(s[-1], s[-2], c='red', s=28)
                fig.savefig(path_to_save_fig + folder + "dapi_spots/smfish_cy5" + f[:-5])
                fig, ax = plt.subplots(1, 1, figsize=(30, 35))
                plt.title(f + " Cy3 green and cy5 red", fontsize=20)
                if img.ndim == 3:
                    ax.imshow(np.amax(img, 0), cmap='gray')
                else:
                    ax.imshow(img, cmap='gray')

                set_cluster_568 = [el[0] for el in positive_cluster_568]
                for s_index in range(len(spots_568)):
                    if labels_568[s_index] in set_cluster_568:
                        s = spots_568[s_index]
                        ax.scatter(s[-1], s[-2], c='green', s=10)

                set_cluster_647 = [el[0] for el in positive_cluster_647]
                for s_index in range(len(spots_647)):
                    if labels_647[s_index] in set_cluster_647:
                        s = spots_647[s_index]
                        ax.scatter(s[-1], s[-2], c='red', s=10)

                fig.savefig(path_to_save_fig + folder + "dapi_spots/green_red" + f[:-5])
                fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                plt.title(f + " Cy3 green and cy5 red", fontsize=20)
                if img.ndim == 3:
                    ax.imshow(np.amax(img, 0), cmap='gray')
                else:
                    ax.imshow(img, cmap='gray')

                for s_index in range(len(spots_568)):
                        s = spots_568[s_index]
                        ax.scatter(s[-1], s[-2], c='green', s=10)

                for s_index in range(len(spots_647)):
                        s = spots_647[s_index]
                        ax.scatter(s[-1], s[-2], c='red', s=10)

                fig.savefig(path_to_save_fig + folder + "dapi_spots/green_red_full" + f[:-5])
                # Cy3 orange
                fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                plt.title(f + " Cy3 channel, orange", fontsize=20)
                if img.ndim == 3:
                    ax.imshow(np.amax(img, 0), cmap='gray')
                else:
                    ax.imshow(img, cmap='gray')

                set_cluster_568 = [el[0] for el in positive_cluster_568]
                for s_index in range(len(spots_568)):
                    if labels_568[s_index] in set_cluster_568:
                        s = spots_568[s_index]
                        ax.scatter(s[-1], s[-2], c='orange', s=10)

                fig.savefig(path_to_save_fig + folder + "dapi_spots/orange_cy3" + f[:-5])

                # cy5 orange
                fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                plt.title(f + " cy5 channel, orange", fontsize=20)
                if img.ndim == 3:
                    ax.imshow(np.amax(img, 0), cmap='gray')
                else:
                    ax.imshow(img, cmap='gray')

                set_cluster_647 = [el[0] for el in positive_cluster_647]
                for s_index in range(len(spots_647)):
                    if labels_647[s_index] in set_cluster_647:
                        s = spots_647[s_index]
                        ax.scatter(s[-1], s[-2], c='orange', s=10)
                fig.savefig(path_to_save_fig + folder + "dapi_spots/orange_cy5" + f[:-5])

                # point cloud
                if not os.path.exists(path_to_save_fig + folder + "convex_hull/"):
                    os.mkdir(path_to_save_fig + folder + "convex_hull/")
                fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                plt.title(f + " Cy3 channel", fontsize=20)
                if img.ndim == 3:
                    ax.imshow(np.amax(img, 0), cmap='gray')
                else:
                    ax.imshow(img, cmap='gray')

                set_cluster_568 = [el[0] for el in positive_cluster_568] + [el[0] for el in negative_cluster_568]
                for s_index in range(len(spots_568)):
                    if labels_568[s_index] in set_cluster_568:
                        s = spots_568[s_index]
                        ax.scatter(s[-1], s[-2], c='orange', s=10)

                from scipy.spatial import ConvexHull
                for c in set_cluster_568:
                    point_cloud = []
                    for s_index in range(len(spots_568)):
                        if labels_568[s_index] == c:
                            point_cloud.append([spots_568[s_index][2], spots_568[s_index][1]])
                    points = np.array(point_cloud)
                    if len(points) < 3:
                        continue
                    hull = ConvexHull(points)
                    for simplex in hull.simplices:
                        ax.plot(points[simplex, 0], points[simplex, 1], 'c')
                        ax.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'o', mec='r', color='none', lw=1,
                                markersize=10)
                fig.savefig(path_to_save_fig + folder + "convex_hull/Cy3_convexhull" + f[:-5])

                fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                plt.title(f + " Cy5 channel", fontsize=20)
                if img.ndim == 3:
                    ax.imshow(np.amax(img, 0), cmap='gray')
                else:
                    ax.imshow(img, cmap='gray')

                set_cluster_647 = [el[0] for el in positive_cluster_647] + [el[0] for el in positive_cluster_647]
                for s_index in range(len(spots_647)):
                    if labels_647[s_index] in set_cluster_647:
                        s = spots_647[s_index]
                        ax.scatter(s[-1], s[-2], c='orange', s=10)

                from scipy.spatial import ConvexHull
                for c in set_cluster_647:
                    point_cloud = []
                    for s_index in range(len(spots_647)):
                        if labels_647[s_index] == c:
                            point_cloud.append([spots_647[s_index][2], spots_647[s_index][1]])
                    points = np.array(point_cloud)
                    hull = ConvexHull(points)
                    for simplex in hull.simplices:
                        ax.plot(points[simplex, 0], points[simplex, 1], 'c')
                        ax.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'o', mec='r', color='none', lw=1,
                                markersize=10)
                fig.savefig(path_to_save_fig + folder + "convex_hull/Cy5_convexhull" + f[:-5])

                fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                plt.title(f + " Cy5 channel", fontsize=20)
                if img.ndim == 3:
                    ax.imshow(np.amax(img, 0), cmap='gray')
                else:
                    ax.imshow(img, cmap='gray')

                set_cluster_647 = [el[0] for el in positive_cluster_647] + [el[0] for el in positive_cluster_647]
                for s_index in range(len(spots_647)):
                    if labels_647[s_index] in set_cluster_647:
                        s = spots_647[s_index]
                        ax.scatter(s[-1], s[-2], c='orange', s=10)

                from scipy.spatial import ConvexHull
                for c in set_cluster_647:
                    point_cloud = []
                    for s_index in range(len(spots_647)):
                        if labels_647[s_index] == c:
                            point_cloud.append([spots_647[s_index][2], spots_647[s_index][1]])
                    points = np.array(point_cloud)
                    hull = ConvexHull(points)
                    for simplex in hull.simplices:
                        ax.plot(points[simplex, 0], points[simplex, 1], 'c')
                        ax.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'o', mec='r', color='none', lw=1,
                                markersize=10)
                fig.savefig(path_to_save_fig + folder + "convex_hull/Cy5_convexhull" + f[:-5])

                if True: # double point_cloud

                    fig, ax = plt.subplots(1, 1, figsize=(30, 20))
                    plt.title(f + " Cy3_green-Cy5-red channel", fontsize=20)
                    if img.ndim == 3:
                        ax.imshow(np.amax(img, 0), cmap='gray')
                    else:
                        ax.imshow(img, cmap='gray')

                    set_cluster_568 = [el[0] for el in positive_cluster_568] + [el[0] for el in negative_cluster_568]
                    for s_index in range(len(spots_568)):
                        if labels_568[s_index] in set_cluster_568:
                            s = spots_568[s_index]
                            ax.scatter(s[-1], s[-2], c='green', s=35)

                    from scipy.spatial import ConvexHull
                    for c in set_cluster_568:
                        point_cloud = []
                        for s_index in range(len(spots_568)):
                            if labels_568[s_index] == c:
                                point_cloud.append([spots_568[s_index][2], spots_568[s_index][1]])
                        points = np.array(point_cloud)
                        hull = ConvexHull(points)
                        for simplex in hull.simplices:
                            ax.plot(points[simplex, 0], points[simplex, 1], 'c', markersize=10,  lw=10)
                            ax.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'o', mec='r',
                                    color='none', lw=30,
                                    markersize=10)

                    set_cluster_647 = [el[0] for el in positive_cluster_647] + [el[0] for el in positive_cluster_647]
                    for s_index in range(len(spots_647)):
                        if labels_647[s_index] in set_cluster_647:
                            s = spots_647[s_index]
                            ax.scatter(s[-1], s[-2], c='red', s=35)

                    from scipy.spatial import ConvexHull
                    for c in set_cluster_647:
                        point_cloud = []
                        for s_index in range(len(spots_647)):
                            if labels_647[s_index] == c:
                                point_cloud.append([spots_647[s_index][2], spots_647[s_index][1]])
                        points = np.array(point_cloud)
                        hull = ConvexHull(points)
                        for simplex in hull.simplices:
                            ax.plot(points[simplex, 0], points[simplex, 1], 'c', markersize=10,  lw=10)
                            ax.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'o', mec='r', color='red',
                                    lw=30,
                                    markersize=10)
                    fig.savefig(path_to_save_fig + folder + "convex_hull/Cy3_Cy5_convexhull" + f[:-5])

                plt.close("all")
                plt.pause(1)
                del fig
                del ax
                gc.collect()


#%%

if __name__ == '__main__':

    print("in the main")
    import torch
    print(f"cuda is available {torch.cuda.is_available()}")
    parser = argparse.ArgumentParser(description='test')

    parser.add_argument('-ptz', "--path_to_czi_folder",
                        type=str,
                        default="/media/tom/250822/czi_folder/",
                        help='path to the parent folder containing the czi')
    parser.add_argument("--list_folder", nargs="+", default=['00_Macrophages/'],
                        help=' list of folders in the czi folders to analyse')
    parser.add_argument('--new_probe', type=str, nargs='+', action='append', default=[],
                        help=" command to add new probes or change parameters of existing one to add it do  "
                        "'--new_probe p1 epsi overlapping  --new_probe p2 20 0.3 ...' "
                        "where  'epsi' is the parameter of the dbscan 'overlapping' "
                        "is the percent of overlap to make a cell positive to a probe")
    parser.add_argument('--manual_threshold_cy3',
                        type=str, default='{"02_NI1230_Lamp3-Cy5_Pdgfra-Cy3_04.tiff": 45}',
                         help=' write a json like the : {"02_NI1230_Lamp3-Cy5_Pdgfra-Cy3_08.tiff": 8,'
                        ' "01_IR5M1236_Lamp3-Cy5_Pdgfra-Cy5_04.tiff": 7} '
                        'to set manually the rna spot detection threshold of possible problematic Cy3 images')
    parser.add_argument('--manual_threshold_cy5', type=str,  default='{"02_NI1230_Lamp3-Cy5_Pdgfra-Cy3_08.tiff": 8,'
                                                                     ' "01_IR5M1236_Lamp3-Cy5_Pdgfra-Cy5_04.tiff": 7}',
                        help="same than manual_threshold_cy3")
    parser.add_argument('-dns', "--dict_name_save", type=str,
                        default="dico_241122",
                        help='key word use to name .npy files storing results of the analysis')

    # cellpose arg default param are for 3D
    parser.add_argument('-d', "--diameter", type=float, default=None, help='cellpose parameter')
    parser.add_argument('-ft', "--flow_threshold", type=float, default=0.75, help='cellpose parameter')
    parser.add_argument('-d3', "--do_3D", type=bool, default=False, help='cellpose parameter')
    parser.add_argument('-m', "--mip", type=bool, default=False, help='cellpose parameter')
    parser.add_argument('-st', "--stitch_threshold", type=float, default=0.4, help='cellpose parameter')
    parser.add_argument('-er', "--erase_solitary", type=int, default=1, help='erase too small nuclei')

    # task to do
    parser.add_argument('-prczi', "--prepare_czi", type=int, default=1, help='if 1 do : prepare_czi to tiff ')
    parser.add_argument('-sg', "--segmentation", type=int, default=1, help='if 1do segmentation ')
    parser.add_argument("--spot_detection", type=int, default=1, help='if 1do spots detection ')
    parser.add_argument("--classify", type=int, default=1, help='if 1do classification / cell type mapping')
    parser.add_argument("--save_plot", type=int, default=1, help='if 1 do plot')

    # parameter for clustering and cell type calling
    parser.add_argument("--scale", nargs='+', default=[300, 103, 103], help='')

    parser.add_argument("--epsi_cluster_cy3", default="é", help="value set in the code if not a float")
    parser.add_argument("--epsi_cluster_cy5", default="e", help="value set in the code if not a float ")
    parser.add_argument("--remove_overlaping", type=int, default=1, help='')
    parser.add_argument("--overlapping_cy3", default="e", help='')
    parser.add_argument("--overlapping_cy5", default="e", help='')
    parser.add_argument("--gpu", type=int, default=1, help='')
    parser.add_argument("--kn_568", type=int, default=3,  help='In pixel size of the xy axis ')
    parser.add_argument("--kn_647", type=int, default=3,  help='In pixel size of the xy axis ')
    parser.add_argument("--port", default=39948)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')



    args = parser.parse_args()

    if args.path_to_czi_folder[-1] != "/":
        args.path_to_czi_folder += "/"

    print(args)
    print()
    print(args.list_folder)
    main(args)

