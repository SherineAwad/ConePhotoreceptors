import scanpy as sc
import anndata as ad
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('obj1')
parser.add_argument('obj2')
parser.add_argument('merged') 
args = parser.parse_args()

obj1 = args.obj1 
obj2 = args.obj2 
merged = args.merged

adata1 = sc.read_h5ad(obj1)
adata2 = sc.read_h5ad(obj2)

merged_adata = ad.concat([adata1, adata2], join="inner")  # "inner" keeps only common genes

print(merged_adata.obs["sample"].unique())

merged_adata.write(merged)

