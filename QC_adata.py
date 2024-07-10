#!/usr/bin/python

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
from utils import tresholds, layer_palette, apply_percentile_treshold, apply_upper_treshold

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

input_path = sys.argv[1]
output_path = sys.argv[2]
slab = sys.argv[3]

# load space ranger output
adata = sc.read_h5ad(input_path)

#adata.var_names_make_unique()
adata.var['mt'] = adata.var.gene_name.str.startswith('MT-')
adata.var['hb'] = adata.var.gene_name.str.contains(("^HB[AB]"))
adata.var['ribo'] = adata.var.gene_name.str.contains(("^RP[LS]"))

# draw plots show raw data quality
p=adata.copy()
sc.pp.calculate_qc_metrics(p, qc_vars=['mt','hb','ribo'],percent_top=None, log1p=False, inplace=True)
sc.pl.violin(p, keys = ['total_counts', 'n_genes_by_counts','pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'], jitter=0.4, rotation= 45, multi_panel=True, size=2.5, save=f'_{slab}_raw_qc.png')

# QC spots
adata=p.copy()
apply_percentile_treshold(adata, 'total_counts')
apply_percentile_treshold(adata, 'n_genes_by_counts', right_bound=False)
apply_upper_treshold(adata, 'pct_counts_mt')
apply_upper_treshold(adata, 'pct_counts_ribo')
apply_upper_treshold(adata, 'pct_counts_hb')
adata.obs['qc_total_counts'].astype(bool).value_counts()
d = {'True': True, 'False': False}
adata.obs['qc_total_counts'].replace(d, inplace=True)
adata.obs['qc_n_genes_by_counts'].replace(d, inplace=True)
adata.obs['qc_pct_counts_mt'].replace(d, inplace=True)
adata.obs['qc_pct_counts_ribo'].replace(d, inplace=True)
adata.obs['qc_pct_counts_hb'].replace(d, inplace=True)
adata.obs['qc_good_spots'] = (adata.obs.qc_total_counts.astype(bool) * adata.obs.qc_n_genes_by_counts.astype(bool) * adata.obs.qc_pct_counts_mt.astype(bool) * adata.obs.qc_pct_counts_ribo.astype(bool) * adata.obs.qc_pct_counts_hb.astype(bool)).astype('string')
vc = adata.obs.groupby('qc_good_spots')['label'].value_counts()
print(f"filter sample number N={adata.obs['qc_good_spots'].value_counts()['False']}\n")
print(vc[vc > 0]['False'])
good_spots = adata.obs[(adata.obs.qc_good_spots == 'True') & (adata.obs.label != 'Empty spots') & (adata.obs.label.notna())].index
adataered = adata[good_spots].copy()
adataered.obs.drop(['qc_total_counts','qc_n_genes_by_counts','qc_pct_counts_mt', 'qc_good_spots', 'qc_pct_counts_hb', 'qc_pct_counts_ribo'], axis=1, inplace=True)
adataered.write_h5ad(output_path + f"{slab}.h5ad")

# draw plots show clean data quality
p=adataered.copy()
sc.pp.calculate_qc_metrics(p, qc_vars=['mt','hb','ribo'],percent_top=None, log1p=False, inplace=True)
sc.pl.violin(p, keys = ['total_counts', 'n_genes_by_counts','pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'], jitter=0.4, rotation= 45, multi_panel=True, size=2.5, save=f'_{slab}_clean_qc.png')
