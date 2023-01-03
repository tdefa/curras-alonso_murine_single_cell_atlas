

# Command line tool to infer cell types using smFISH data with two channels of two different marker genes.

- This code is specifically designed to analyse the images from *A murine single cell atlas of the lung response to radiation injury*, A.Curras-Alonso et al.

## System requirements

- Tested on **Python 3.8** running on **Ubuntu 20.04** on a **Dell XPS-15**.
- Images are stored in the czi format.

## Installation guide

1. We recommend using miniconda from [**here**](https://docs.conda.io/en/latest/miniconda.html)
2. Create your conda environment and all required dependencies with the provided file `ENV.txt` (takes usually a few minutes): `conda list --explicit>ENV.txt`

## Demo

demo.sh

## Instructions for use

### Runing the cell type mapping algorithm

1. Execute the following command in the code directory
2. Cell type mapping plots will be store in the "Figure" folder.
3. Replace the parameters by your configuration.  
4. run  `python main_cluster.py --help` for more information
5. Results are saved in a .npy file dictionary

**Example**:

- Image names must contain the name of the gene that was imaged.
- Code contains specific analysis parameters for each gene, e.g. the minimal distance for the dbscan.
- For genes that don't contain hard-coded parameters an option allows to provide them as additional input (see below)

``` bash
python main_cluster.py --path_to_czi_folder /media/tom/250822/czi_folder/  
--list_folder experiment1/ experiment4/
--new_probe Pdgfratest 35 0.42
--new_probe Pdgfratest2 40 0.62
--dict_name_save analysis2022  
```

**Parameters**:
- `--path_to_czi_folder`: path to the parental folder containing subfolders with czi images of different experiments
- `--list_folder`: names of folders in the parental containing experiments to be anaylzed
- `--new_probe`: permits to specify analysis parameters for genes that are not listed in the code: `espi` parameter in dbscan (minimal distance between points in a cluster), minimal overlapping between marker gene and nucleus to make a nucleus positive to this marker
- `dict_name_save`: key word use to name .npy files storing results of the analysis

## Generating Excel files with analysis results

1. Execute the following commands in the code directory indicating the name of the .npy file dictionary where analysis results are saved <br />
- The script main_generate_excel_one_probe.py will generate an exels file for each gene, containing for each image, the number of cell positive to this genes and the aproximative mean cytoplasm and nuclei size of positive cell . 

- The script main_generate_excel_two_probe.py will generate an exels file for each pair of genes, containing for each image, the number of cell positive to the pair of genes and the aproximative mean cytoplasm and nuclei size of positive cell. 

**Example**:
```bash
python  main_generate_excel_one_probe.py \ 
--path_save /media/tom/Elements1/test_pipeline/exels/ \
--path_folders_czi /media/tom/250822/czi_folder/   \
--list_folder experience1/ experience4/ \
--list_probes Lamp3 \ 
--list_probes 'Cap' 'aCap' 'CEC' 'acap' \ 
--dict_name_save analysis2022.npy \ 
```

**Parameters**:
- `--path_to_czi_folder`: Path to the folder where the exels containing cell type calling results will be saved.
- `--path_folders_czi`: Path to the parental folder containing subfolders with czi images of different experiments.
- `--list_folder`: Names of folders in the parental containing experiments to be anaylzed.
- `--list_probes`: Probes to analyse, add space to separate probes that have many names.
- `--dict_name_sav`: Name of the .npy file dictionary storing cell type calling results from main_cluster.py




**Example**:
```bash
python  main_generate_excel_two_probe.py \ 
--path_save /media/tom/Elements1/test_pipeline/exels/ \ 
--path_folders_czi /media/tom/250822/czi_folder/ 
--list_folder experience1/ experience4/ 
--list_probes1 Lamp3 \ #[First set of probes to analyse] <br />
--list_probes2 'Cap' 'aCap' 'CEC' 'acap' \ #[second set of probes that will be compared the the first set, add space to separate probes that have many names (like here)]  <br />
--dico_stat_name analysis2022.npy \ 


```
**Parameters**:Same than for  main_generate_excel_one_probe except than  `--list_probes` is replaced by:

- `--list_probes1`: Name(s) of the first probe to analyse
- `--list_probes2`: Name(s) of the second probe to analyse 

Then an excel will be generated with cell positive to both (probe1, probe2)


#### Documentation

use python [name of the command] --help for a more extensive documentation

