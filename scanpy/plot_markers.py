import scanpy as sc
import sys
import importlib_metadata
import anndata 
import argparse
import scanpy.external as sce
import pandas as pd

sys.modules['importlib.metadata'] = importlib_metadata


parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject

combined_adata = sc.read(myObject)
parts = myObject.split("_")

fname = parts[1].split(".") 
filename = fname[0] +"markerGenes" +".png" 


marker_genes  = {
        "MG": ["Rlbp1a","Gfap","Apoe", "Hes1", "Notch1", "Aqp4", "Pax6", "Prdx6", "Slc1a3", "Abca8a", "Vim", "Aldh1a", "Hes5"],
        "RPC": ["Sox2", "Ascl1", "Gli1", "Sfrp2", "Fgf15"], 
        "Rod": ["Rho", "Nrl", "Otx2", "Crx", "Guca1b", "Rom1", "Nr2e3"], 
        "Cones": ["Opn1mw", "Opn1sw", "Arr3", "Thrb", "Otx2", "Gnat2", "Cgna3"], 
        "BC": ["Sebox", "Bhlhe23", "Cabp5", "Vsx2", "Prkca", "Pcp4", "Isl1", "Vsx1", "Grm6", "Trpm1"], 
        "AC": ["Gad1", "Gad2", "Slc6a9", "Tfap2b", "Prox1", "Pax6", "Calb2", "Pcp4", "Elavl3", "Elavl4"], 
        "HC": ["Lhx1", "Cbln4", "Calb1", "Nefl", "Nefm", "Onecut1", "Onecut2"], 
        "RGC":["Sncg", "Thy1", "Ebf3", "Rbfox3", "Isl12", "Pou4f1", "Pou4f2", "Pou4f3", "Rbpms"], 
        "Microglia": ["Ptprc", "Cx3cr1", "Csf2rb", "Sall1"], 
        "Percicytes": ["Kcnj8"], 
        "Astrocytes": ["Pax2", "Igf2", "Gfap", "Pdgfra"], 
        "Endothelial": ["Tie1", "Pecam1"],
        "Olignocytes": ["Mbpa"] 
    }

# Filter marker genes to include only those present in the dataset
filtered_marker_genes = {
    cluster: [gene for gene in genes if gene in combined_adata.var_names]
    for cluster, genes in marker_genes.items()
}

# Remove empty clusters (where no genes were found in the dataset)
filtered_marker_genes = {k: v for k, v in filtered_marker_genes.items() if v}

# Plot dotplot only if there are valid genes
if filtered_marker_genes:
    sc.pl.dotplot(combined_adata, filtered_marker_genes, groupby="leiden", standard_scale="var", save=filename)

# Plot UMAP scatter plots for each gene in the filtered dictionary
for cluster, genes in filtered_marker_genes.items():
    for gene in genes:
        sc.pl.scatter(combined_adata, color=gene, title=f'{cluster} - {gene}', basis='umap', save=f'_{gene}.png')

















