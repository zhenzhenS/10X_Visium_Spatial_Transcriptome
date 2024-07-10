def norm_by_layer(adata):
  for i in np.unique(adata.obs.sample_id.values):
    mean_value=[]
    mean_value_df=[]
    for k in np.unique(adata[adata.obs.sample_id == i].obs.layer.values):
      mean_value.append(adata[(adata.obs.sample_id == i) & (adata.obs.layer == k)].X.mean(axis=0))
    mean_value_df=pd.DataFrame(np.concatenate(mean_value),columns=adata[(adata.obs.sample_id == i) & (adata.obs.layer == k)].var_names) #all layer mean
    all_mean=mean_value_df.mean(axis=0).to_numpy() #all layer mean calculate sample mean
    adata[adata.obs.sample_id == i].X=csr_matrix(adata[adata.obs.sample_id == i].X-all_mean)
  return adata
