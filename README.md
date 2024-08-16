This repository contains the script to automatically annotate plaque scRNAseq datasets with out atlas. 

How to use it

Download neccesary files:

python download.py

Install neccesary packages:

Python:
- scanpy
- scarches
- anndata
- numpy
- pandas
- matplotlib
- scipy
- rpy2
- anndata2ri

R:
- scran

Alternativly, we also provide a docker image with all the packages installed and files downloaded at: https://github.com/matmu/cell_type_annotation. 

We also provide a fully automated solution in an interactive web framework at: https://www.archmap.bio/#/


Run script:

python automatic_mapping.py <filename> <lognorm> <celltype> <ensembleIDs> 
where the <filename> should be the path to the h5ad object that you want to map, <lognorm> to specify if your counts are log-normalized (true/false), <celltype> to specify if your object has a "cell_type_level2" column filled with "unknown" (true/false) and <ensemblIDs> to specify if your varnames are ensemblIDs (true/false). We recommend to input the counts after cells are removes with QC and to use ensemblIDs if possible. It you only have gene names, the script tries to map the gene names to ensemblIDs.

Example:
We also provide a example file that we also used to validate in the publication.

python automatic_mapping.py example_data/Bashore_postQC_noCITE_noR.h5ad false false false

Because its quite a lot of cells (75k cells) the normalization step will take some time. You will also need a Nvidia GPU (cuda) to run the script.