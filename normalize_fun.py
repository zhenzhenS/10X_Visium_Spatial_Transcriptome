def norm_spots(adata):
  for i in np.unique(adata.obs.sample_id.values):
    mean_value=[]
    mean_value_df=[]
    for k in np.unique(adata[adata.obs.sample_id == i].obs.layer.values):
      mean_value.append(adata[(adata.obs.sample_id == i) & (adata.obs.layer == k)].X.mean(axis=0))
    mean_value_df=pd.DataFrame(np.concatenate(mean_value),columns=adata[(adata.obs.sample_id == i) & (adata.obs.layer == k)].var_names) #all layer mean
    all_mean=mean_value_df.mean(axis=0).to_numpy() #all layer mean calculate sample mean
    adata[adata.obs.sample_id == i].X=csr_matrix(adata[adata.obs.sample_id == i].X-all_mean)
  return adata

def norm_layer(adata):
    # calculate mean gene expr for each sample
    gene_mean_list = []
    sample_layer_list = adata.obs.sample_id.unique().tolist()
    for sample in sample_layer_list:
        gene_mean_list.append(adata[adata.obs.sample_id == sample].X.mean(axis=0).reshape(-1, 1))
    # convert to DataFrame
    gene_mean_df = pd.DataFrame(np.concatenate(gene_mean_list, axis=1), columns=sample_layer_list, index=adata.var_names)
    gene_mean_df.head()
    
    # subtract mean and return adata
    for sample in sample_layer_list:
        columns = adata.obs[adata.obs.sample_id == sample].index.to_list()
        for column in columns:
            adata[column].X = adata[column].X - gene_mean_df.loc[:, sample].values
    return adata
