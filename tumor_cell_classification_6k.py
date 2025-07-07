import os
import numpy as np
import pandas as pd
import gseapy as gp
from sklearn.semi_supervised import LabelSpreading


def softmax(x):
    e_x = np.exp(x - np.max(x)) 
    return e_x / e_x.sum()


def calc_simp_score(df, qq = .66):
    data = df.to_numpy().astype(np.float32)
    softmax_data = np.apply_along_axis(softmax, axis = 1, arr = data)
    
    row_max = np.max(softmax_data, axis = 1)
    row_mean_ex_max = (softmax_data.sum(axis = 1) - row_max) / (softmax_data.shape[1] - 1)
    
    new_col = row_max - row_mean_ex_max
    df_transformed = pd.DataFrame(softmax_data, columns = df.columns, index = df.index)
    df_transformed['simp_score'] = new_col
    
    th = df_transformed['simp_score'].quantile(qq)
    df_transformed['is_top_qq'] = df_transformed['simp_score'] >= th
    df_transformed['argmax_label'] = np.argmax(df_transformed.iloc[:,:4].values, axis = 1)
    
    mapping = dict(zip(np.arange(4), df.columns))
    df_transformed['argmax_st'] = df_transformed['argmax_label'].map(mapping)
    
    return df_transformed



def classify_tumor_cells(sample, cts, sigs, ks, n_cpu):
    
    gsea_res_by_fov = {}
    unique_fovs = cts.fov.unique()

    for ff in unique_fovs:
        tmp = cts[cts['fov'] == ff].iloc[:,:-1]
        ss_res = gp.ssgsea(
            tmp.T,
            gene_sets = sigs,
            threads = n_cpu,
            random_seed = 37,
            verbose = True,
            min_size = 0
        )
        print(f'done fov {ff}')
        gsea_res_by_fov[ff] = ss_res
        
    to_concat = []
    for val in gsea_res_by_fov.values():
        tmp = val.res2d.copy()
        tmp = tmp.pivot(index = 'Name', columns = 'Term', values = 'NES')
        tmp = calc_simp_score(tmp)
        to_concat.append(tmp)

    res = pd.concat(to_concat)
    res['partial_label'] = res['argmax_label'].values
    res.loc[~res['is_top_qq'], 'partial_label'] = -1
    
    for k in ks:
        ls = LabelSpreading(kernel = 'knn', n_neighbors = k, n_jobs = n_cpu)
        ls.fit(X = cts.iloc[:,:-1].values, y = res['partial_label'].values)
        res[f'ls_{k}'] = ls.transduction_
        print(f'finished ls k = {k}')
        
    return res