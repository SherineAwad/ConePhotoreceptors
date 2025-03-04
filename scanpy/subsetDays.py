import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import scipy.sparse as sp
import numpy as np



sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')

args = parser.parse_args()

myObject =  args.myObject

combined_adata = sc.read_h5ad(myObject, backed="r")


combined_adata.obs['day_group'] = combined_adata.obs['sample'].apply(lambda x: '15days' if '15day' in x else ('35days' if '35day' in x else 'Other'))


adata_15days = combined_adata[combined_adata.obs['day_group'] == '15days']

adata_35days = combined_adata[combined_adata.obs['day_group'] == '35days']

adata_15days = adata_15days.to_memory()
adata_35days = adata_35days.to_memory()

adata_15days = adata_15days.copy()
adata_35days = adata_35days.copy()

adata_15days.write("15days.h5ad", compression="gzip")
adata_35days.write("35days.h5ad", compression="gzip")
