#!/usr/bin/python
##################################################
#   read visium data, add label and sample_id,   #
#   save as h5ad format file.                    #
##################################################

# import packages
import igraph
import torch
import sklearn
import scanpy as sc
import anndata as ann
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

input_path = sys.argv[1]
ann_path = sys.argv[2]
slab = sys.argv[3]
output_path = sys.argv[4]

# load space ranger output
adata = sc.read_visium(input_path)
adata.var['gene_name']=adata.var.index
adata.var.index=adata.var.gene_ids

# load layer assignment
ann = pd.read_csv(ann_path, index_col=0)
ann['sample_id']=slab 
ann.columns=['label','sample_id']

# add layer and other label into adata
adata.obs['label']=ann.label
adata.obs['sample_id']=ann.sample_id

adata = adata[adata.var.gene_name.notna()]
adata.obs['layer']=adata.obs.label
adata.obs.layer.replace({'L1R': 'L1', 'L1L': 'L1','L2R': 'L2', 'L2L': 'L2','L3R': 'L3', 'L3L': 'L3',
  'L4R': 'L4', 'L4L': 'L4','L5R': 'L5', 'L5L': 'L5','WMR': 'WM', 'WML': 'WM','L6aR': 'L6', 'L6bR': 'L6',
  'L6aL': 'L6', 'L6bL': 'L6','L6a':'L6','L6b':'L6','6a':'L6'}, inplace=True)
print(adata)
print(adata.obs)
print(adata.var)
print(adata.obs.layer.value_counts())
adata.write_h5ad(output_path + f'raw/{slab}_raw.h5ad')
