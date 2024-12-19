# Plaque scRNAseq Dataset Annotation

This repository contains a script to automatically annotate plaque scRNAseq datasets using our atlas.

## How to Use

### Install Necessary Packages

#### Python:
- gdown
- scanpy
- scarches
- anndata
- numpy
- pandas
- matplotlib
- scipy
- rpy2
- anndata2ri

#### R:
- scran

### Download Necessary Files

Run the following command:

```
python download.py
```

Alternatively, we provide a Docker image with all packages installed and files downloaded at: [https://github.com/matmu/cell_type_annotation](https://github.com/matmu/cell_type_annotation).

We also offer a fully automated solution in an interactive web framework at: [https://www.archmap.bio/#/](https://www.archmap.bio/#/)

### Run Script

Use the following command to use the level 1 annotations (recommended):

```
python automatic_mapping_level1.py <filename> <lognorm> <celltype> <ensembleIDs>
```

or level 2 annotations:

```
python automatic_mapping_level2.py <filename> <lognorm> <celltype> <ensembleIDs>
```

Where:
- `<filename>`: Path to the h5ad object you want to map
- `<lognorm>`: Specify if your counts are log-normalized (true/false)
- `<celltype>`: Specify if your object has a "cell_type_level2" column filled with "unknown" (true/false)
- `<ensembleIDs>`: Specify if your varnames are ensemblIDs (true/false)

We recommend inputting the counts after cells are removed with QC and using ensemblIDs if possible. If you only have gene names, the script attempts to map the gene names to ensemblIDs using a mapping used in our atlas.

### Example

We provide an example file that we used for validation in the publication:

```
python automatic_mapping_level1.py example_data/Bashore_postQC_noCITE_noR.h5ad false false false
```

Note: Because it contains a large number of cells (75k), the normalization step will take some time. You will also need an Nvidia GPU (CUDA) to run the script.
