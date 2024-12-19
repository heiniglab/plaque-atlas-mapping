#!/usr/bin/env python

import scanpy as sc
from scarches.models.scpoli import scPoli
import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import issparse
import scipy.sparse as sparse
from scarches.models.base._utils import _validate_var_names
import sys
import logging
import os

################## input arguments ##################

def str_to_bool(value):
    return value.lower() in ('true', '1', 't', 'yes', 'y')

adata_input_file = sys.argv[1]
lognorm_bool = str_to_bool(sys.argv[2])
cell_type_bool = str_to_bool(sys.argv[3])
ensembl_bool = str_to_bool(sys.argv[4])
#tissue_name_file = sys.argv[5]

# remove the .h5ad ending
#tissue_name = tissue_name_file.split(".")[0]
# remove the "TS_" prefix
#tissue_name = tissue_name[3:]

#print(f"Create dir healthy_mapping/TS/{tissue_name}2")
#os.makedirs(f"healthy_mapping/TS/{tissue_name}2")


################## Loading in the atlas and model ##################

print("Load model and atlas....")

adata = sc.read_h5ad("data/Big-Atlas-level122-log1p-hvg.h5ad")

#load trained reference model
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
    
scpoli_model = scPoli(
adata=adata,
condition_keys="sample",
cell_type_keys="cell_type_level2", 
embedding_dims=10,
recon_loss='mse',
)

scpoli_loaded = scpoli_model.load(dir_path="models/reference_retraining22", adata=adata)


################## Loading in the mapping data ##################

print("Loading the mapping data...")
#adata_bashore = sc.read_h5ad("tmp/Bashore_postQC_noCITE_noR.h5ad")
adata_bashore = sc.read_h5ad(adata_input_file)

# get correct layer for Tabula
#adata_bashore.X = adata_bashore.layers['decontXcounts']
#adata_bashore.obs.rename(columns={'donor': 'sample'}, inplace=True)

#adata_bashore.X = adata_bashore.X.toarray()
adata_bashore_full = adata_bashore.copy()

adata_bashore.obs["sample"] = [sample + "_query" for sample in adata_bashore.obs["sample"]]

#map gene ids to ensembl
if ensembl_bool == False:

    print("Mapping and aggregation of ensembl genes...")
    ensembl_id_df = pd.read_csv("data/gene_names_to_ensembl_ALLFOUND_allfernandez_no6_withallslysz.csv")
    gene_to_ensembl = dict(zip(ensembl_id_df['gene_name'], ensembl_id_df['ensembl_id']))
    # Map the variable names in AnnData
    adata_bashore.var['original_gene_names'] = adata_bashore.var_names
    adata_bashore.var_names = [gene_to_ensembl[gene] if gene in gene_to_ensembl else gene for gene in adata_bashore.var_names]

    non_ENSG_vars = adata_bashore.var_names[~adata_bashore.var_names.str.startswith('ENSG')]
    # remove not mapped genes
    # Convert non_ENSG_vars to a set for faster lookup
    non_ENSG_vars_set = set(non_ENSG_vars)

    # Filter out the variables that are in non_ENSG_vars_set
    adata_bashore = adata_bashore[:, ~adata_bashore.var_names.isin(non_ENSG_vars_set)]


    #aggregation
    # Convert the sparse matrix (if it is sparse) to a dense DataFrame
    adata_df = pd.DataFrame(adata_bashore.X.toarray() if issparse(adata_bashore.X) else adata_bashore.X, 
                            index=adata_bashore.obs_names, 
                            columns=adata_bashore.var_names)

    # Group by gene names and sum the counts
    aggregated_data = adata_df.groupby(adata_df.columns, axis=1).sum()

    # Prepare the new 'var' DataFrame, keeping the first occurrence of each gene
    unique_var = adata_bashore.var.loc[~adata_bashore.var.index.duplicated(keep='first')]

    # Create a new AnnData object with aggregated data
    adata_agg = anndata.AnnData(X=aggregated_data, obs=adata_bashore.obs, var=unique_var.loc[aggregated_data.columns])

    # 'adata_agg' now has unique gene names and aggregated counts

    adata_agg.X = sparse.csr_matrix(adata_agg.X)


    #checks
    print("Check how many cells have zero counts for all genes...")
    cellwise_sum = adata_agg.X.sum(axis=1)
    num_cells_zero_counts = (cellwise_sum == 0).sum()
        
    if num_cells_zero_counts>0:
        print(num_cells_zero_counts, " cells were found with 0 counts across all genes! Removing these cells now...")
        adata_agg = adata_agg[cellwise_sum > 0, :]
    
    adata_final = adata_agg.copy()

else:
    adata_final = adata_bashore.copy()


if lognorm_bool == False:
    print("Log normalize data...")
    import rpy2.rinterface_lib.callbacks
    from rpy2.robjects import pandas2ri
    import anndata2ri
    import rpy2.robjects as ro

    # Copy the data
    adata_pp = adata_final.copy()

    # Normalize data
    print("Normalize the data...")
    sc.pp.normalize_total(adata_pp, target_sum=1e6)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, svd_solver="arpack")
    sc.pp.neighbors(adata_pp, n_pcs=30)
    sc.tl.leiden(adata_pp, key_added='groups', resolution=0.22)

    # Setup rpy2 and suppress R warnings
    rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

    # Automatically convert rpy2 outputs to pandas dataframes
    pandas2ri.activate()
    anndata2ri.activate()

    # Preprocess variables for scran normalization
    input_groups = adata_pp.obs['groups']
    data_mat = adata_final.X.T.toarray()

    # Define the R script as a string
    r_script = """
    library(scran)
    size_factors = calculateSumFactors(data_mat, clusters=input_groups, min.mean=0.1)
    """

    # Execute the R script
    ro.globalenv['data_mat'] = data_mat
    ro.globalenv['input_groups'] = input_groups
    ro.r(r_script)
    size_factors = ro.globalenv['size_factors']

    # Remove intermediate adata_pp
    del adata_pp

    # Normalize adata
    adata_final.obs['size_factors'] = size_factors
    adata_final.X /= adata_final.obs['size_factors'].values[:, None]
    sc.pp.log1p(adata_final)


if cell_type_bool == False:
    adata_final.obs['cell_type_level2'] = "unknown"


print("Select the varnames for the model....")
varnames_path = "models/reference_retraining22/var_names.csv"
var_names = np.genfromtxt(varnames_path, delimiter=",", dtype=str)
adata2 = _validate_var_names(adata_final, var_names)


################## Run the model ##################

### INSERT YOUR .H5AD ADATA HERE. has to have a "sample" obs and ensembl id vars, and normalized + log1p, cell_type_level2 obs with "unknown" aswell. 
### also varnames need to be selected
adata_query2 = adata2.copy() 

scpoli_query2 = scPoli.load_query_data(
    adata=adata_query2,
    reference_model=scpoli_loaded,
    labeled_indices=[],
)


scpoli_query2.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10
)


#scpoli_query2.save("models/bashore-mapping", overwrite=True)
#scpoli_query2 = scpoli_query2.load(dir_path="models/bashore-mapping", adata=adata_query2)

adata_query2.X = adata_query2.X.astype("float32")

results_dict = scpoli_query2.classify(adata_query2, scale_uncertainties=True)


#get whole atlas including reference
#get latent representation of reference data
scpoli_query2.model.eval()
data_latent_source = scpoli_query2.get_latent(
    adata,
    mean=True
)

adata_latent_source = sc.AnnData(data_latent_source)
adata_latent_source.obs = adata.obs.copy()

#get latent representation of query data
data_latent= scpoli_query2.get_latent(
    adata_query2,
    mean=True
)

adata_latent = sc.AnnData(data_latent)
adata_latent.obs = adata_query2.obs.copy()

#get label annotations
adata_latent.obs['cell_type_pred'] = results_dict['cell_type_level2']['preds'].tolist()
adata_latent.obs['cell_type_uncert'] = results_dict['cell_type_level2']['uncert'].tolist()
#adata_latent.obs['classifier_outcome'] = (adata_latent.obs['cell_type_pred'] == adata_latent.obs['cell_type_level1'])

#join adatas
adata_latent_full = adata_latent_source.concatenate(
    [adata_latent],
    batch_key='query'
)


#adata_latent_full.write("healthy_mapping/Hu/Hu_test1.h5ad")

adata_latent_full.obs['cell_type_pred_ref'] = np.where(
    adata_latent_full.obs['query'].isin(['0']),  
    "Reference",                                
    adata_latent_full.obs['cell_type_pred']     
)

#adata_latent_full.write("healthy_mapping/Hu/Hu_test2.h5ad")



################## Postprocessing ##################

print("Start post processing")
adata_latent_full.obs['cell_type_level2_all'] = [
    row['cell_type_level2'] if row['cell_type_pred_ref'] == 'Reference' else row['cell_type_pred'] 
    for _, row in adata_latent_full.obs.iterrows()
]

# rename cell_type_level2_all to cell_type_level2 and remove cell_type_level2_all
del adata_latent_full.obs['cell_type_level2']
adata_latent_full.obs['cell_type_level2'] = adata_latent_full.obs['cell_type_level2_all']
del adata_latent_full.obs['cell_type_level2_all']

bashore = adata_latent_full[adata_latent_full.obs['cell_type_pred_ref'] != "Reference"]

#delete leftover obs
del bashore.obs["conditions_combined"]
del bashore.obs["cell_type_pred"]
del bashore.obs["cell_type_pred_ref"]
del bashore.obs["query"]
del bashore.obs["n_counts"]
del bashore.obs["dataset"]
del bashore.obs["cell_type_level1"]

#remove the added suffix in barcodes to have the same input and output barcodes
bashore.obs.index = [idx[:-2] for idx in bashore.obs.index]

# add cell type preds and uncertainties to original adata object

# Extract the observations from bashore
bashore_obs = bashore.obs[['cell_type_level2', 'cell_type_uncert']]

# Create a DataFrame for adata_bashore_full with the same indices
adata_bashore_full_obs = pd.DataFrame(index=adata_bashore_full.obs_names)

# Merge the data from bashore_obs into adata_bashore_full_obs
# If the barcode is not found in bashore, the merged result will have NaN
merged_obs = adata_bashore_full_obs.merge(bashore_obs, left_index=True, right_index=True, how='left')

# Fill NaN values with "unknown"
merged_obs.fillna('unknown', inplace=True)

# Assign the merged observations back to adata_bashore_full
adata_bashore_full.obs[['cell_type_level2', 'cell_type_uncert']] = merged_obs


#save output file with embeddings
print("Saving h5ad files...")
sc.pp.neighbors(bashore, n_neighbors=15)
sc.tl.umap(bashore)

bashore.write("output/embedding_level2.h5ad")
adata_bashore_full.write("output/full_level2.h5ad")


############ PLOTTING ##############


color_palette_level2 = {
    'B cell': '#2ca02c',                       # green

    'CD4 T cell': '#6495ED',                   # Cornflower Blue
    'CD8 T cell': '#0047AB',                   # Cobalt Blue

    'EndoMT EC': '#FFD700',                    # Gold
    'Pro-Angiogenic EC': '#ED9121',            # Carrot Orange

    'Fibroblast': '#e377c2',                   # bright pink
    'Fibromyocyte': '#f7b6d2',                 # pastel pink
    'Smooth Muscle Cell': '#7b4173',           # deeper purple

    'Other Macrophage': '#ece2d0',              # very light brown
    'TREM2+/Foamy Macrophage': '#c4a484',       # light tan-brown
    'HMOX1+ Macrophage': '#bf9b7a',             # medium brown
    'Inflammatory Macrophage': '#8b4513',       # saddle brown
    'PLIN2+/TREM1+ Macrophage': '#5b3a29',      # dark brown

    'Mast cell': '#d62728',                    # red

    'Monocyte': '#c7c7c7',                     # lighter gray
    'Neutrophil': '#17becf',                   # teal

    'NK cell': '#98df8a',                      # light green

    'Plasma cell': '#9467bd',                  # muted purple

    'Conventional dendritic cell 1': '#d4b9da',  # light violet
    'Conventional dendritic cell 2': '#807dba',  # medium violet
    'Plasmacytoid dendritic cell': '#4a1486'     # deep violet
}




# Plot UMAP with custom color palette
sc.settings.figdir = "output/" # if you want to change the output path for figures
sc.pl.umap(
    bashore,
    color='cell_type_level2',
    palette=color_palette_level2,
    show=True,
    frameon=False,
    title = "Level 2",
    #save=f"{tissue_name}_level1.png"
    save=f"level2.png"
)

sc.pl.umap(
    bashore,
    color='cell_type_uncert',
    show=True,
    frameon=False,
    title = "Uncertainties",
    #save=f"{tissue_name}_uncertainties.png"
    save=f"uncertainties_level2.png"
)

