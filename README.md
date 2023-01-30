

# Command line tool to infer cell types using smFISH data with two channels of two different marker genes.

- This code is specifically designed to analyse the images from the study

    *A murine single cell atlas of the lung response to radiation injury*, S.Curras-Alonso et al.

## System requirements

- Tested on **Python 3.8** running on **Ubuntu 20.04** on a **Dell XPS-15**.
- Images are stored in the czi format.

## Installation guide

1. We recommend using miniconda from [**here**](https://docs.conda.io/en/latest/miniconda.html).
2. Create your conda environment and all required dependencies with the provided file `ENV.yml` (takes usually a few minutes): `conda env create -n ENVNAME --file ENV.yml`

## Demo

To run the analysis on a small demo data set: 

1. Download demo data from [here](https://doi.org/10.5281/zenodo.7501932).
2. Run the provided bash file in this repository in terminal with enabled conda environment: `./demo.sh` (modify `path` in `demo.sh` if needed)

## Instructions for use

### Runing the cell type mapping algorithm

1. Execute the following command in the code directory.
2. Cell type mapping plots will be store in the "Figure" folder.
3. Replace the parameters by your configuration.  
4. run  `python main_cluster.py --help` for more information.
5. Results are saved in a .npy file dictionary.

**Example**:

- Image names must contain the name of the gene that was imaged.
- Code contains specific analysis parameters for each gene, e.g. the minimal distance for the dbscan.
- For genes that don't contain hard-coded parameters an option allows to provide them as additional input (see below)

``` bash
python main_cluster.py --path_to_czi_folder /media/tom/250822/czi_folder/ \
--list_folder experiment1/ experiment4/ \
--new_probe Pdgfratest 35 0.42 \
--new_probe Pdgfratest2 40 0.62 \
--dict_name_save analysis2022 \
```

**Parameters**:
- `--path_to_czi_folder`: path to the parental folder containing subfolders with czi images of different experiments
- `--list_folder`: names of folders in the parental containing experiments to be anaylzed
- `--new_probe`: permits to specify analysis parameters for genes that are not listed in the code: name if the probe, `espi` parameter in dbscan (minimal distance between points in a cluster), minimal overlapping between marker gene and nucleus to make a nucleus positive to this marker. <br /> 
Genes with pre-defined default values : 'Lamp3', 'Pdgfra', 'Chil3', 'Cap', 'aCap', 'acap', 'Ptprb', 'Ptprb1', 'Fibin', 'C3ar1', 'Hhip', 'Mki67', 'Serpine1', 'Apln', 'Pecam1', 'CEC', 'Rtkn2'
- `dict_name_save`: key word use to name .npy files storing results of the analysis

## Generating Excel files with analysis results

Commands have to be executed in the code directory, options are explained in more detail below.


### Result files for one gene

The script `main_generate_excel_one_probe.py` will generate an Excel file for the specified gene, 
containing for each image, the number of positive cells for this gene and the estimated mean 
cellular and nuclear volume. 

**Example**:
```bash
python  main_generate_excel_one_probe.py \ 
--path_save /media/tom/Elements1/test_pipeline/excel/ \
--path_folders_czi /media/tom/250822/czi_folder/   \
--list_folder experience1/ experience4/ \
--list_probes 'Lamp3' \ 
--list_probes 'Cap' 'aCap' 'CEC' 'acap' \ 
--dict_name_save analysis2022.npy \ 
```

**Parameters**:
- `--path_to_czi_folder`: Path to the folder where the exels containing cell type calling results will be saved.
- `--path_folders_czi`: Path to the parental folder containing subfolders with czi images of different experiments.
- `--list_folder`: Names of folders in the parental containing experiments to be anaylzed.
- `--list_probes`: Probes to analyse, add space to separate probes that have many names.
- `--dict_name_sav`: Name of the .npy file dictionary storing cell type calling results from main_cluster.py

### Result files for pairs of genes

The script `main_generate_excel_two_probes.py` will generate an Excel file for all specified pairs of genes, 
containing for each image, the number positive cells for each gene, and estimated mean 
cellular and nuclear volume. . 

**Example**:
```bash
python  main_generate_excel_two_probe.py \
--path_save /media/tom/Elements1/test_pipeline/exels/ \
--path_folders_czi /media/tom/250822/czi_folder/ \
--list_folder experience1/ experience4/ \
--list_probes1 'Lamp3' \
--list_probes2 'Cap' 'aCap' 'CEC' 'acap' \
--dico_stat_name analysis2022.npy \

```
**Parameters**:Same than for  main_generate_excel_one_probe except than  `--list_probes` is replaced by:

- `--list_probes1`: Name(s) of the first gene to analyse
- `--list_probes2`: Name(s) of the second genes to analyse 

Then an excel will be generated with cell positive to both (probe1, probe2)


#### Documentation

Use `python [name of the command] --help` for a more extensive documentation.

