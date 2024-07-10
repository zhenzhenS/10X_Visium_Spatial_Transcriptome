#!/usr/bin/python

import sklearn
import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys

input_path = sys.argv[1]
slab = sys.argv[2]
output_path = sys.argv[3]

obj_list=[]
for i in os.listdir(input_path):
        if slab in i:
                tmp=sc.read_h5ad(input_path + i)
                obj_list.append(tmp)
adata=sc.concat(obj_list,merge="same",uns_merge='unique')
adata.obs_names_make_unique()
adata.write_h5ad(output_path + f"{slab}.h5ad")
print(adata)
print(adata.obs)
print(adata.var)
print(adata.uns['spatial'].keys())
