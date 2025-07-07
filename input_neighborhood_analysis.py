import os
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, coo_matrix
import pyarrow.feather as pf


def read_adj(path_to_files, smp, radius = 500):
    mat = pf.read_feather(os.path.join(path_to_files, f'{smp}_adj_mat_{str(radius)}px.feather'))
    dimnames = pf.read_feather(os.path.join(path_to_files, f'{smp}_adj_dimnames.feather'))['cell_id'].values
    row_idx = mat['row'].values
    col_idx = mat['col'].values
    values = mat['val'].values
    n = row_idx.max() + 1
    sp_mat = coo_matrix((values, (row_idx, col_idx)), shape = (n, n)).tocsr()
    return sp_mat



def calc_neighborhood(mdata, obs_key, adj_mat, drop_categories):
    
    labs = mdata[obs_key].dropna().unique().tolist()
    labs = [l for l in labs if l not in drop_categories]
    target_labs = [l for l in labs if l.startswith(('GPM', 'MTC', 'NEU', 'PPR'))]
    source_labs = [l for l in labs if l not in target_labs]

    cell_ids = mdata.index.to_numpy()
    cell_id_to_idx = {cid: i for i, cid in enumerate(cell_ids)}
    labels = mdata.loc[cell_ids, obs_key].values

    source_mask = np.isin(labels, source_labs)
    target_mask = np.isin(labels, target_labs)

    source_idxs = np.where(source_mask)[0]
    target_idxs = np.where(target_mask)[0]
    target_labels = labels[target_mask]

    target_label_mat = np.zeros((len(target_idxs), len(target_labs)), dtype = np.uint8) #1hot encoding of the subtypes
    
    label_to_col = {lab: i for i, lab in enumerate(target_labs)}
    
    for i, lab in enumerate(target_labels):
        target_label_mat[i, label_to_col[lab]] = 1

    submat = sp_mat[source_idxs, :][:, target_idxs]  # shape: (num_sources, num_targets)

    #(num_sources × num_targets) @ (num_targets × num_labels)
    res = submat @ target_label_mat  # shape: (num_sources, num_labels)
    res = pd.DataFrame(true_res, index = cell_ids[source_mask], columns = [f'{st}' for st in target_labs])
    
    return res