# Plaque scRNAseq Dataset Annotation

This repository contains a script to automatically annotate plaque scRNAseq datasets using our atlas (preprint: https://www.biorxiv.org/content/10.1101/2024.09.11.612431v1).

## How to Use

### Install Necessary Packages

#### Python:
- gdown (5.2.0)
- scanpy (1.9.6)
- scarches (0.5.10)
- anndata (0.10.3)
- numpy (1.23.4)
- pandas (2.2.1)
- matplotlib (3.8.2)
- scipy (1.11.4)
- rpy2 (3.5.14)
- anndata2ri (1.3.1)

#### R:
- scran (1.30.0)

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
python automatic_mapping_level1.py example_data/Hu_subset.h5ad false false true
```

Note: Because it contains a large number of cells (37k), the normalization step will take some time. You will also need an Nvidia GPU (CUDA) to run the script.
