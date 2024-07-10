#!/usr/bin/python

import sklearn
import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys

input_path = sys.argv[1]
adata=sc.read_h5ad(input_path)
adata=adata[adata.obs.sample_id !="chimp_13302"]
adata.obs["log_lib_size"] = np.log(adata.obs["total_counts"])
sc.tl.pca(adata)
sc.pl.pca(adata,
          color=['layer', 'sample_id', 'log_lib_size'],
          size=200,ncols=2, show=True,save=f'pca_{sys.argv[2]}.pdf')
