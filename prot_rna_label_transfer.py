import os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from sklearn.preprocessing import MinMaxScaler
import anndata as ad
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment


def project_coords_in_common_space(rna_coords, prot_coords, scale = True):
    mtx1 = rna_coords.copy()
    mtx2 = prot_coords.copy()

    mtx1 -= np.mean(mtx1, 0)
    mtx2 -= np.mean(mtx2, 0)
    
    if scale:
        scaler = MinMaxScaler()

        mtx1 = scaler.fit_transform(mtx1)
        mtx2 = scaler.fit_transform(mtx2)
    
    return mtx1, mtx2


def match_by_corr_and_dist(rna_adata, prot_adata, fov_mapping, prot_rna_mapping, spatial_key, classification_key, verbose = False):
    res = []
    for rna_fov, prot_fov in fov_mapping.items():
        if ',' in rna_fov:
            rna_fov = rna_fov.split(',')
        else:
            rna_fov = [rna_fov]
        rna_tmp = rna_adata[rna_adata.obs.fov.isin(rna_fov)].copy()
        rna_counts = rna_tmp.X.toarray()
        rna_coords = pd.DataFrame(rna_tmp.obsm[spatial_key])
        if isinstance(prot_fov, str):
            prot_fov = [prot_fov]
        prot_tmp = prot_adata[prot_adata.obs.fov.isin(prot_fov)].copy()
        prot_counts = prot_tmp.X.toarray()
        prot_coords = pd.DataFrame(prot_tmp.obsm[spatial_key])
        prot_coords[1] = -prot_coords[1] #conserve correct orientation of spatial coords in fov
        
        rna_coords_transformed, prot_coords_transformed = project_coords_in_common_space(rna_coords, prot_coords)
        
        correlations = np.nan_to_num(1 - cdist(rna_counts, prot_counts, metric = 'correlation'), nan = -1)
        distances = cdist(rna_coords_transformed, prot_coords_transformed)
        gamma = -np.log(1e-3)/np.max(distances)
        
        util = correlations * (1 + np.exp(-gamma * distances))
        
        row_ind, col_ind = linear_sum_assignment(util, maximize = True)
        
        if verbose:
            print(f'optimization completed for rna fovs {rna_fov} and prot fovs {prot_fov}')
        df = pd.DataFrame({
            'rna_id' : rna_tmp.obs_names[row_ind],
            'prot_id' : prot_tmp.obs_names[col_ind],
            'classification' : rna_tmp[rna_tmp.obs_names[row_ind]].obs[classification_key].values,
            'util' : util[row_ind, col_ind]
        })
        res.append(df)
    return pd.concat(res, axis = 0)