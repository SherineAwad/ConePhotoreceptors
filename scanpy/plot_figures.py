import scanpy as sc
import sys
import importlib_metadata
import anndata 
import argparse
import scanpy.external as sce
import numpy as np
import matplotlib.pyplot as plt
sys.modules['importlib.metadata'] = importlib_metadata


parser = argparse.ArgumentParser()
parser.add_argument('myObject')
parser.add_argument('genesFile') 
args = parser.parse_args()

myObject =  args.myObject
genesFile = args.genesFile 

combined_adata = sc.read(myObject)
parts = myObject.split(".")
sample_condition = parts[0]

figureN = sample_condition.split(".") 

filename = parts[0] +"_" + figureN[0] +".png" 


def read_genes_from_file(file_path):
    with open(file_path, 'r') as file:
        genes = [line.strip() for line in file.readlines()]
    return genes

#genes_of_interest = read_genes_from_file(genesFile)

#def diff_expression_dotplot(adata, genes_of_interest, sample_condition):
#    adata.obs['condition'] = adata.obs['sample'].apply(lambda x: 'Ctrl' if 'S1' in x else 'KO')
#    sc.tl.rank_genes_groups(adata, groupby='condition', method='t-test')
#    de_results_dict = adata.uns['rank_genes_groups']
#    de_results_df = pd.DataFrame({
#        'genes': de_results_dict['names'][group] for group in de_results_dict['names'].dtype.names
#    })
#    sc.pl.dotplot(adata, var_names=genes_of_interest, groupby='celltype', save=filename )
#diff_expression_dotplot(combined_adata, genes_of_interest, sample_condition)

'''
#One gene in a plot 
combined_adata.obs["condition"] = combined_adata.obs["sample"].apply(lambda x: 'Ctrl' if 'S1' in x else 'KO')
combined_adata = combined_adata[~combined_adata.obs["celltype"].isna(), :]
combined_adata = combined_adata[~combined_adata.obs["condition"].isna(), :]
combined_adata.obs["celltype"] = combined_adata.obs["celltype"].astype(str)
combined_adata.obs["condition"] = combined_adata.obs["condition"].astype(str)
filename= sample_condition + "_Gls.png"
sc.pl.dotplot(combined_adata, var_names=["Gls"], groupby="celltype", 
              standard_scale="var", swap_axes=True, 
              figsize=(10, 8), save=filename)

'''





combined_adata.obs['condition'] = combined_adata.obs['celltype'].apply(lambda x: 'KO' if 'KO' in x else 'Ctrl')

sc.pl.dotplot(combined_adata, var_names=["Gls"], groupby="celltype", 
              hue="condition", standard_scale="var", swap_axes=True, 
              save="Gls_dotplot.png")




