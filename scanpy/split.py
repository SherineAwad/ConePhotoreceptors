import scanpy as sc
import sys
import importlib_metadata
import anndata 
import argparse
import scanpy.external as sce

sys.modules['importlib.metadata'] = importlib_metadata


parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject

combined_adata = sc.read(myObject)

combined_adata.var.head()
combined_adata.obs.head()

combined_adata.obs['sample'] = [
    f"15dayS{barcode.split('-')[1]}" if barcode.split('-')[1] in ['1', '2'] else "Other"
    for barcode in combined_adata.obs.index
]


combined_adata.write(myObject)







