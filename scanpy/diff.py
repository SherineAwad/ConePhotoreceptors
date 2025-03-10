import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import scipy
import scipy.sparse as sp
import numpy as np
import pandas as pd
import h5py


sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject

parts = myObject.split(".")
newObject = "diff_" + myObject

sample = parts[0] 
combined_adata = sc.read_h5ad(myObject, backed="r")

combined_adata = combined_adata.to_memory()
print("Cell Metadata (obs):")
print(combined_adata.obs.head())  # First few rows of cell metadata

print("\nGene Metadata (var):")
print(combined_adata.var.head())  # First few rows of gene metadata

#sc.pp.normalize_total(combined_adata, target_sum=1e4)

sc.pp.log1p(combined_adata)
combined_adata.X = np.nan_to_num(combined_adata.X, nan=0, posinf=0, neginf=0)
#print("Replaced NaN and Inf values with 0.")


print("Shape of combined_adata.X:", combined_adata.X.shape)
print("Number of NaNs in X:", np.isnan(combined_adata.X).sum())
print("Number of Infs in X:", np.isinf(combined_adata.X).sum())
print("Min value in X:", np.min(combined_adata.X))
print("Max value in X:", np.max(combined_adata.X))

non_zero_genes = (combined_adata.X.sum(axis=0) > 0)  # No `.A1` needed
combined_adata = combined_adata[:, non_zero_genes]
print(f"Filtered genes with all-zero values. New shape: {combined_adata.shape}")



sc.pp.highly_variable_genes(combined_adata)

sc.tl.rank_genes_groups(combined_adata, groupby="sample", method="wilcoxon")




de_results_dict = combined_adata.uns['rank_genes_groups']

de_results_df = pd.DataFrame({
    "group": [],
    "gene": [],
    "logfoldchange": [],
    "pval": [],
    "pval_adj": []
})

filename = sample + "dge.csv" 
for group in de_results_dict['names'].dtype.names:
    for i, gene in enumerate(de_results_dict['names'][group]):
        de_results_df = de_results_df.append({
            "group": group,
            "gene": gene,
            "logfoldchange": de_results_dict['logfoldchanges'][group][i],
            "pval": de_results_dict['pvals'][group][i],
            "pval_adj": de_results_dict['pvals_adj'][group][i]
        }, ignore_index=True)

de_results_df.to_csv(filename, index=False)

combined_adata.obs["Pathway_Score"] = np.random.rand(combined_adata.n_obs)
sc.pl.umap(combined_adata, color="Pathway_Score",save="pathway.png")

combined_adata.write(newObject, compression="gzip")

