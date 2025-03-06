import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
parts = myObject.split("_")  
newObject = "clustered_" + parts[1] 

combined_adata = sc.read(myObject)
fname = parts[1].split(".") [0]  

figure_name = fname+ "_clusters.png"
sc.tl.pca(combined_adata)
sc.pp.neighbors(combined_adata)
sc.tl.leiden(combined_adata, n_iterations=3)  #to adjust resolution add ,resolution=1.5)
sc.tl.umap(combined_adata)
sc.pl.umap(combined_adata, color=["leiden"], save= figure_name,legend_loc="on data")


combined_adata.write(newObject,compression="gzip")
