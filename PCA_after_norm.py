#!/usr/bin/python

import sklearn
import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys

input_path = sys.argv[1]
slab = sys.argv[2]

adata=sc.read_h5ad(input_path)
print(adata.obs.label.value_counts())
sc.pl.violin(adata,
             keys = ['n_genes_by_counts', 'total_counts','pct_counts_ribo', 'pct_counts_mt'],
             jitter=0.4, groupby = 'sample_id', rotation= 45,save=f'{slab}_raw_sample_id.pdf')
sc.pl.violin(adata,
             keys = ['n_genes_by_counts', 'total_counts','pct_counts_ribo', 'pct_counts_mt'],
             jitter=0.4, groupby = 'label', rotation= 90,save=f'{slab}_raw_label.pdf')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata)
sc.pl.pca(adata,
          color=['label', 'sample_id'],
          ncols=2, show=True,save=f'_{slab}_raw.pdf')
