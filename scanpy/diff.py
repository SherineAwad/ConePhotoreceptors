import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import scipy.sparse as sp
import numpy as np
import pandas as pd



sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject

parts = myObject.split(".")
newObject = "diff_" + parts[0]

sample = parts[0] 

combined_adata = sc.read_h5ad(myObject, backed="r")

combined_adata = combined_adata.to_memory()
sc.pp.normalize_total(combined_adata, target_sum=1e4) 
sc.pp.log1p(combined_adata)
sc.pp.highly_variable_genes(combined_adata)
sc.tl.rank_genes_groups(combined_adata, groupby="sample", method="wilcoxon")
print(combined_adata.uns["rank_genes_groups"].keys())
print(combined_adata.uns["rank_genes_groups"]["names"])


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



combined_adata.write(newObject, compression="gzip")
