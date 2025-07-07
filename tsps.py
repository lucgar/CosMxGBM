import os
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad
import scanpy as sc
import squidpy as sq



def calc_proximity_score(adata, classification_key, categories, spatial_key, radius = 250, score_key = 'tsps', inplace = True):
    
    tmp_adata = adata[adata.obs[classification_key].isin(categories)]
    
    adj, _ = sq.gr.spatial_neighbors(
    tmp_adata,
    spatial_key = spatial_key,
    coord_type = 'generic',
    radius = radius,
    copy = True
    )
    
    labels = tmp_adata.obs[classification_key].values
    
    row_indices, col_indices = adj.nonzero()
    same_label_mask = labels[row_indices] == labels[col_indices]
    
    same_label_mat = csr_matrix((same_label_mask.astype(int), (row_indices, col_indices)), shape = adj.shape)
    
    num = same_label_mat.sum(axis = 1).A1
    den = adj.sum(axis=1).A1
    
    tsps = np.divide(num, den, out = np.full_like(num, np.nan, dtype = float), where = den != 0)
    tsps = pd.Series(tsps, index = tmp_adata.obs_names)
    
    adata.obs[score_key] = tsps
    
    if not inplace:
        return tsps
    
    
    
def calc_proximity_score_with_permut(
    adata,
    classification_key,
    categories,
    spatial_key,
    radius = 250,
    n_perms = 100,
    random_state = 37
):

    tmp_adata = adata[adata.obs[classification_key].isin(categories)]

    adj, _ = sq.gr.spatial_neighbors(
        tmp_adata,
        spatial_key=spatial_key,
        coord_type='generic',
        radius=radius,
        copy=True
    )

    labels = tmp_adata.obs[classification_key].values
    obs_names = tmp_adata.obs_names

    row_indices, col_indices = adj.nonzero()
    
    # Calc stat
    same_label_mask = labels[row_indices] == labels[col_indices]
    same_label_mat = csr_matrix((same_label_mask.astype(int), (row_indices, col_indices)), shape=adj.shape)
    num = same_label_mat.sum(axis = 1).A1
    den = adj.sum(axis = 1).A1
    tsps = np.divide(num, den, out = np.full_like(num, np.nan, dtype = float), where=den != 0)
    
    results = pd.DataFrame(index = obs_names)
    results['tsps'] = tsps

    # Permut
    for i in range(n_perms):
        np.random.seed(random_state * i)
        shuffled_labels = np.random.permutation(labels)
        boot_mask = shuffled_labels[row_indices] == shuffled_labels[col_indices]
        boot_mat = csr_matrix((boot_mask.astype(int), (row_indices, col_indices)), shape=adj.shape)
        num_boot = boot_mat.sum(axis = 1).A1
        tsps_boot = np.divide(num_boot, den, out=np.full_like(num_boot, np.nan, dtype = float), where= den != 0)
        results[f'tsps_boot{i+1}'] = tsps_boot

    return results