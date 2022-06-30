
#%% utils for sort data to plot

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
    print(key_cell_name)
    print(gene_smfish)
    raise Exception("Probe not detected: please the probe name should be in the fiel name")


def get_mouse_name(image_name):

    if any(word in image_name for word in ["NI", 'Ctrl']):
        if all(word in image_name for word in ["1225"]):
            return "1225"
        elif all(word in image_name for word in ["1230"]):
            return "1230"
        else:
            return "2323"

    if any(word in image_name for word in ['IR5M']):
        if "1236" in image_name:
            return "1236"
        elif "1249" in image_name:
            return "1249"
        elif "1250" in image_name:
            return "1250"
        elif "2330" in image_name:
            return "2330"
        else:
            return "2201"
