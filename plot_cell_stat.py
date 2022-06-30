import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from scipy.stats import mannwhitneyu

#%%



#%%
def get_list_NI_IRM_by_animal(dataframe_per_files,
                              column_numerateur = "nb_positive",
                              NI_name =  ["1225", '1230', '2323'],
                              ir5m_name = ['1236', '1249', '1250', '2330', "2201"],
                              gy10_name  = ["1251", "1259", "1260"]):


    if column_numerateur == "nb_positive":
        dataframe_per_files = dataframe_per_files[["experiment", "image_name", "nb_nuclei", column_numerateur]]
        dataframe_per_files = dataframe_per_files[dataframe_per_files.nb_nuclei != ""]
        dataframe_per_files = dataframe_per_files.dropna()
        dataframe_per_files['tagret'] = dataframe_per_files[column_numerateur].astype(int) / dataframe_per_files[
            "nb_nuclei"].astype(int) * 100


    if column_numerateur == "average_point_cloud_size" or column_numerateur ==  "average_nuclei_size":
        dataframe_per_files = dataframe_per_files[dataframe_per_files.average_point_cloud_size != "not defined"]
        dataframe_per_files = dataframe_per_files[dataframe_per_files.average_nuclei_size != "not defined"]
        dataframe_per_files = dataframe_per_files[dataframe_per_files.average_point_cloud_size != "Not define"]
        dataframe_per_files = dataframe_per_files[dataframe_per_files.average_nuclei_size != "Not define"]
        dataframe_per_files = dataframe_per_files[["experiment", "image_name", "nb_nuclei", "nb_positive", column_numerateur]]
        dataframe_per_files = dataframe_per_files[dataframe_per_files.nb_nuclei != 0]
        dataframe_per_files = dataframe_per_files[dataframe_per_files.nb_positive > 0]
        dataframe_per_files = dataframe_per_files.dropna()

        dataframe_per_files['tagret'] = dataframe_per_files[column_numerateur].astype(int)


    if column_numerateur == "nb_positive_both":
        dataframe_per_files = dataframe_per_files[["experiment", "image_name", "nb_nuclei", column_numerateur]]
        dataframe_per_files = dataframe_per_files[dataframe_per_files.nb_nuclei != ""]
        dataframe_per_files = dataframe_per_files.dropna()
        dataframe_per_files['tagret'] = dataframe_per_files[column_numerateur].astype(int)



    dico_qm = {}

    dico_qm['NI'] = {}
    for name in NI_name:
        if name == '2323':
            print( "waaaarrring here check mouse name")
            #print(dataframe_per_files['image_name'])
            #print(name)
            NI_name = dataframe_per_files.loc[
            (dataframe_per_files['image_name'].str.contains('NI') | dataframe_per_files['image_name'].str.contains('Ctrl'))
            &
            -dataframe_per_files['image_name'].str.contains('1230')
            & -dataframe_per_files['image_name'].str.contains('1225'), ['tagret']]
        else:
            NI_name = dataframe_per_files.loc[
            (dataframe_per_files['image_name'].str.contains('NI') | dataframe_per_files['image_name'].str.contains('Ctrl'))
            &
            dataframe_per_files['image_name'].str.contains(name), ['tagret']]

        dico_qm['NI'][name] = list(NI_name['tagret'])
        print( NI_name['tagret'])

    dico_qm['IR5M'] = {}
    for name in  ir5m_name:
        if name == "2201":
            print( "waaaarrring here check mouse name")
            #print(dataframe_per_files['image_name'])
            #print(name)

            IRM_name = dataframe_per_files.loc[dataframe_per_files['image_name'].str.contains('IR5M_') &
                                               -dataframe_per_files['image_name'].str.contains('1236')
                                               & -dataframe_per_files['image_name'].str.contains('1249')
                                               & -dataframe_per_files['image_name'].str.contains('1250')
                                               & -dataframe_per_files['image_name'].str.contains('2330'), [
                                                   'tagret']]
            if len(IRM_name['tagret']) == 0:
                IRM_name = dataframe_per_files.loc[dataframe_per_files['image_name'].str.contains('IR5M2201') &
                                                   dataframe_per_files['image_name'].str.contains(name), ['tagret']]

        else:
            IRM_name = dataframe_per_files.loc[dataframe_per_files['image_name'].str.contains('_IR5M') &
                                               dataframe_per_files['image_name'].str.contains(name), ['tagret']]



        dico_qm['IR5M'][name] = IRM_name['tagret']
    dico_qm['IR5M10gy'] = {}
    for name in gy10_name:
        print(name)
        IRM10gy_name = dataframe_per_files.loc[dataframe_per_files['image_name'].str.contains('IR5M10Gy') &
                                           dataframe_per_files['image_name'].str.contains(name), ['tagret']]
        dico_qm['IR5M10gy'][name] = IRM10gy_name['tagret']


    #IRM = dataframe_per_files.loc[
    #    dataframe_per_files['image_name'].str.contains('IR1M|IR2M|IR3M|IR4M'), ['tagret']]


    return dico_qm

#%%
def plot_scatter(dico_qm, dico_color, dico_label,
                 title="average size of point cloud per sample for ",
                 axis_y_label="size of point cloud in  μm3",
                 path_to_save = ""):
    """
    list_NI is a list of list [list_value, name]
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    np.random.seed(seed=40)
    x_tick = []
    x_tick_name = []
    x_labels = []
    x_index = 1
    xmax_NI = 0
    xmax_IR5M10gy = 0
    if len(np.concatenate(list(dico_qm["NI"].values()))) > 0:
        for mouse in  dico_qm["NI"]:
            if len( dico_qm["NI"][mouse]) == 0:
                continue
            y =  dico_qm["NI"][mouse]
            x = [x_index] * len( dico_qm["NI"][mouse]) + np.random.rand(len( dico_qm["NI"][mouse])) - 0.5
            s, c = np.array([2.0] * len( dico_qm["NI"][mouse])), np.array([dico_color[mouse]] * len( dico_qm["NI"][mouse]))
            s *= 10.
    
            ax.scatter(x, y, s, c, label=dico_label[mouse])
            x_tick.append(x_index)
            x_labels.append("NI " + mouse)

    if len(np.concatenate(list(dico_qm["NI"].values()))) > 0:
        list_m = np.concatenate([dico_qm["NI"][k] for k in dico_qm["NI"]])
        median_NI = np.median(list_m)
        mean_NI = np.mean(list_m)
        xmax_NI = x_tick[-1]
        ax.hlines(y=median_NI, xmin=0.5, xmax=xmax_NI + 0.5, linewidth=3)
        print("NI MEDIAN %s" % median_NI)
        x_index += 1.3
        x_tick_name.append("NI")

    if len(np.concatenate(list(dico_qm["IR5M10gy"].values())))  > 0:

        for mouse in dico_qm["IR5M10gy"]:
            if len(dico_qm["IR5M10gy"][mouse]) == 0:
                continue
            y = dico_qm["IR5M10gy"][mouse]
            x = [x_index] * len(dico_qm["IR5M10gy"][mouse]) + np.random.rand(len(dico_qm["IR5M10gy"][mouse])) - 0.5
            s, c = np.array([2.0] * len(dico_qm["IR5M10gy"][mouse])), np.array(
                [dico_color[mouse]] * len(dico_qm["IR5M10gy"][mouse]))
            s *= 10.

            ax.scatter(x, y, s, c, label=dico_label[mouse])
            x_tick.append(x_index)
            x_labels.append("IR5M10gy " + mouse)

    if  len(np.concatenate(list(dico_qm["IR5M10gy"].values())))  > 0:
        list_m = np.concatenate([dico_qm["IR5M10gy"][k] for k in dico_qm["IR5M10gy"]])
        median_IR5M10gy = np.median(list_m)
        mean_IR5M10gy = np.mean(list_m)
        xmax_IR5M10gy = x_tick[-1]
        ax.hlines(y=median_IR5M10gy, xmin=max(0.5,xmax_NI + 0.5), xmax=xmax_IR5M10gy + 0.5, linewidth=3)
        print("IR5M10gy MEDIAN %s" % median_IR5M10gy)
        x_index += 1.3
        x_tick_name.append("IR5M10gy")

    xmax_IR5M10gy = max(xmax_IR5M10gy, xmax_NI)
    if len(np.concatenate(list(dico_qm["IR5M"].values()))) > 0:

        for mouse in dico_qm["IR5M"]:
            if len(dico_qm["IR5M"][mouse]) == 0:
                continue
            y = dico_qm["IR5M"][mouse]
            x = [x_index] * len(dico_qm["IR5M"][mouse]) + np.random.rand(len(dico_qm["IR5M"][mouse])) - 0.5
            s, c = np.array([2.0] * len(dico_qm["IR5M"][mouse])), np.array(
                [dico_color[mouse]] * len(dico_qm["IR5M"][mouse]))
            s *= 10.

            ax.scatter(x, y, s, c, label=dico_label[mouse])
            x_tick.append(x_index)
            x_labels.append("IR5M " + mouse)

    if  len(np.concatenate(list(dico_qm["IR5M"].values()))) > 0:
        list_m = np.concatenate([dico_qm["IR5M"][k] for k in dico_qm["IR5M"]])
        median_IR5M = np.median(list_m)
        mean_IR5M = np.mean(list_m)
        xmax_IR5M = x_tick[-1]
        ax.hlines(y=median_IR5M, xmin=max(0.5,xmax_IR5M10gy + 0.5), xmax=xmax_IR5M + 0.5, linewidth=3)
        print("IR5M MEDIAN %s" % median_IR5M)
        x_index += 1.3
        x_tick_name.append("IR5M17gy")







    # draw vertical line from (70,100) to (70, 250)
    ax.set_xticks([i+1 for i in range(len([k  for k in dico_qm.keys() if len(np.concatenate(list(dico_qm[k].values()))) > 0]))])

    ax.set_xticklabels(x_tick_name, rotation=45)
    ax.tick_params(axis='x', which='major', labelsize=20)
    ax.tick_params(axis='y', which='major', labelsize=15)
    ax.set_ylabel(axis_y_label, fontsize=15)
    fig.suptitle(title, fontsize=20)
    try:
        list_NI = np.concatenate(list(dico_qm["NI"].values()))
        list_10gy = np.concatenate(list(dico_qm["IR5M10gy"].values()))
        list_IRM = np.concatenate(list(dico_qm["IR5M"].values()))

        U1, pn17 = mannwhitneyu(list_NI,
                             list_IRM)

        U1, pn10 = mannwhitneyu(list_NI,
                             list_10gy)

        U1, pn1017 = mannwhitneyu(list_10gy,
                                list_IRM)


        ax.set_title('P value  of the Mann-Whitney test NI against IR5M (17gy): ' + str(round(pn17, 6)).ljust(6, '0')
                      + '\n' + 'P value  NI against IR5M (10gy): ' + str(round(pn10, 6)).ljust(6, '0')
                     +  '\n' + 'P value   IR5M (17gy) against IR5M (10gy): ' + str(round(pn1017, 6)).ljust(6, '0') ,
                     fontsize=15)
    except Exception as e:
        ax.set_title(str(e),
                     fontsize=15)
    ax.legend(bbox_to_anchor=(1.05, 1),
              loc='upper left')
    plt.show()

    p=None

    fig.savefig(path_to_save, bbox_inches='tight')

    return #round(mean_NI, 4), round(median_NI, 4), round(mean_IR5M, 4), round(median_IR5M, 4), p


#%%


if __name__ == "__main__":

    #%%

    dico_color = {
        "1225": '#1f77b4',
        "1230": '#ff7f0e',
        "2323": '#2ca02c',

        "1251": '#207a93',
        "1259": '#d645d5',
        "1260": '#1b0d39',

        "1236": '#d62728',
        "1249": '#9467bd',
        "1250": '#8c564b',
        "2201": '#e377c2',
        "2330": '#7f7f7f',
    }

    dico_label = {
        "1225": "NI_1225",
        "1230": 'NI_1230',
        "2323": 'NI_2323',

        "1251": '#IR5M_10gy_1251',
        "1259": '#IR5M_10gy_1259',
        "1260": '#IR5M_10gy_1260',

        "1236": 'IR5M_17gy_1236',
        "1249": 'IR5M_17gy_1249',
        "1250": 'IR5M_17gy_1250',
        "2201": 'IR5M_17gy_2201',
        "2330": 'IR5M_17gy_2330',
    }



    ##one cell
    path_to_exels = "/home/tom/Bureau/sandra050122/exels_all/one_cells/"
    path_to_save = "/home/tom/Bureau/sandra050122/plot_all/nuclei_vol/"

    list_probes = [ ['Lamp3'],  ['Pdgfra'] ,['Serpine1'], ['Ptprb'], ["Rtkn2"],
                 ['Apln'], ['Chil3'],  ['Fibin'], ['C3ar1'],
                 ['Hhip'], ['Mki67'], ['Pecam1'],  ['Cap', 'aCap', 'CEC'],]
  #  list_probes =  [ ['Chil3']]
    dico_cellname_probe = {
        'Lamp3': "AT2",
        'Pdgfra': "fibroblast",
        "Serpine1": "scenescent",
        'Ptprb': "vessel EC",
        'Apln': "capillary EC with Alpn",
        'Chil3': "AM",
        "Fibin": "capillary EC with Fibin",
        'C3ar1': "IM",
        'Hhip': "myifibroblast",
        'Mki67': "cycling cells",
        'Pecam1': "EC",
        'Cap': "capillary EC",
        'aCap': "capillary EC",
        'CEC': "capillary EC",
        "Rtkn2": "AT1"
    }

    for prb in list_probes:

        try:
            dataframe_per_files = pd.read_excel(path_to_exels  + prb[0]+ ".xls")
        except:
            continue

        if np.sum([type(i) == type(" ") for i in list(dataframe_per_files["image_name"])]) == 0:
            continue




        dico_qm = get_list_NI_IRM_by_animal(dataframe_per_files,
                                NI_name=["1225", '1230', '2323'],
                                  ir5m_name=['1236', '1249', '1250', '2330', "2201"],
                                  gy10_name=["1251", "1259", "1260"],
                                column_numerateur ='average_nuclei_size'# "average_point_cloud_size" #,
                                )

        if len([k  for k in dico_qm.keys() if len(np.concatenate(list(dico_qm[k].values()))) > 0]) == 0:
            print(prb[0])
            continue


        plot_scatter(dico_qm, dico_color, dico_label,
                         title= "average nuclei volume per sample %s (probe %s +)" % (dico_cellname_probe[prb[0]], prb[0]),
                         axis_y_label="average nuclei volume in µm3",
                         path_to_save = path_to_save  +dico_cellname_probe[prb[0]]+ "_"+ prb[0])


#%%
    ##one cell
    path_to_exels = "/home/tom/Bureau/sandra050122/exels_all/one_cells/"
    path_to_save = "/home/tom/Bureau/sandra050122/plot_all/point_cloud_vol/"

    list_probes = [ ['Lamp3'],  ['Pdgfra'] ,['Serpine1'], ['Ptprb'], ["Rtkn2"],
                 ['Apln'], ['Chil3'],  ['Fibin'], ['C3ar1'],
                 ['Hhip'], ['Mki67'], ['Pecam1'],  ['Cap', 'aCap', 'CEC'],]
  #  list_probes =  [ ['Chil3']]
    dico_cellname_probe = {
        'Lamp3': "AT2",
        'Pdgfra': "fibroblast",
        "Serpine1": "scenescent",
        'Ptprb': "vessel EC",
        'Apln': "capillary EC with Alpn",
        'Chil3': "AM",
        "Fibin": "capillary EC with Fibin",
        'C3ar1': "IM",
        'Hhip': "myifibroblast",
        'Mki67': "cycling cells",
        'Pecam1': "EC",
        'Cap': "capillary EC",
        'aCap': "capillary EC",
        'CEC': "capillary EC",
        "Rtkn2": "AT1"}

    for prb in list_probes:

        try:
            dataframe_per_files = pd.read_excel(path_to_exels  + prb[0]+ ".xls")
        except:
            continue

        if np.sum([type(i) == type(" ") for i in list(dataframe_per_files["image_name"])]) == 0:
            continue




        dico_qm = get_list_NI_IRM_by_animal(dataframe_per_files,
                                NI_name=["1225", '1230', '2323'],
                                  ir5m_name=['1236', '1249', '1250', '2330', "2201"],
                                  gy10_name=["1251", "1259", "1260"],
                                column_numerateur ="average_point_cloud_size" #'average_nuclei_size'# "average_point_cloud_size" #,
                                )

        if len([k  for k in dico_qm.keys() if len(np.concatenate(list(dico_qm[k].values()))) > 0]) == 0:
            print(prb[0])
            continue


        plot_scatter(dico_qm, dico_color, dico_label,
                         title= "average point cloud volume per sample %s (probe %s +)" % (dico_cellname_probe[prb[0]], prb[0]),
                         axis_y_label="average point cloud volume in µm3",
                         path_to_save = path_to_save  +dico_cellname_probe[prb[0]]+ "_"+ prb[0])




#%% ### two cells
    from pathlib import Path
    import re
    path_to_save = "/media/tom/Transcend/sandra050122/sandra03022022/plot_all/two_porbes_ppositive_proportion/"
    path_to_exels = "/media/tom/Transcend/sandra050122/sandra03022022/exels_all/pair_cells/"
    path_g = Path("/media/tom/Transcend/sandra050122/sandra03022022/exels_all/pair_cells/Pecam1_Ptprb.xls") # Path(path_to_exels)
    for xxxx in path_g.glob("*.xls"):
        xxxx = Path("/media/tom/Transcend/sandra050122/sandra03022022/exels_all/pair_cells/Pecam1_Apln.xls")
        dataframe_per_files = pd.read_excel(str(xxxx))
        probe1, probe2, _ = str(xxxx).replace('.', '_').split('/')[-1].split('_')
        if np.sum([type(i) == type(" ") for i in list(dataframe_per_files["image_name"])]) == 0:
            print(probe1)
            print(probe2)
            print()
            #continue


        dico_qm = get_list_NI_IRM_by_animal(dataframe_per_files, #column='percent_cell',
                                column_numerateur="nb_positive_both",
                                  NI_name=["1225", '1230', '2323'],
                                 # ir5m_name=['1236', '1249', '1250', '2330', "2201"],
                                            ir5m_name = ['1250'],
                                  gy10_name=["1251", "1259", "1260"])

        if len([k for k in dico_qm.keys() if len(np.concatenate(list(dico_qm[k].values()))) > 0]) == 0:
            print(probe1)
            print(probe2)
            print()
            #continue

        plot_scatter(dico_qm, dico_color, dico_label,
                         title= "percent of of (%s + & %s + )" % (probe1, probe2),
                         axis_y_label="proportion percent",
                         path_to_save = path_to_save  +probe1 + "_" + probe2)
