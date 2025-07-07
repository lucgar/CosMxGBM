import os
import numpy as np
import pandas as pd
import gseapy as gp
from sklearn.semi_supervised import LabelSpreading


def run_ssgsea(gex, signature_dict, random_state = 37):
    #transpose the matrix so that you have genes on rows
    tmp = gex.T.copy()
    tmp_gsea = gp.ssgsea(
    data = tmp,
    gene_sets = signature_dict,
    min_size = 0,
    verbose = True,
    threads = -1,
    seed = random_state
    )
    rr = tmp_gsea.res2d
    ss = rr.pivot(index = 'Term', columns = 'Name', values = 'NES')
    return ss


def simp_score(col):
    if np.all(col <= 0):
        return np.nan
    else:
        col = col.astype('float').copy()
        num = np.exp(col)
        probs = num / np.sum(num)
        sorted_probs = np.sort(probs)
        maxx = sorted_probs[-1]
        mean = np.mean(sorted_probs[:-1])
        return maxx - mean
    

def get_partial_labels(ss, gex, percentile = .66):
    scores = ss.apply(simp_score, axis = 0)
    th = scores.quantile(percentile)
    confident_classif = ss.loc[:,scores[scores >= th].index].apply(np.argmax)
    partial_labels = pd.Series(-1, index = gex.index) #### modify the step of the index assignment
    partial_labels[confident_classif.index] = confident_classif
    return partial_labels, scores


def get_cells_by_fov(cts):
    fov_id = pd.Series([int(x.split('_')[2]) for x in cts.index], index = cts.index)
    cells_by_fov = {f'fov_{str(i)}' : cells.tolist() for i, cells in fov_id.groupby(fov_id).groups.items()}
    return cells_by_fov


def ensemble_ls(gex, partial_labels, ks, clamping_factor):
    if (len(partial_labels) <= max(ks)) or np.all(partial_labels == -1):
        return partial_labels
    gex = gex.copy()
    res = np.empty((gex.shape[0], len(ks)))
    for i, k in enumerate(ks):
        ls = LabelSpreading(
        kernel = 'knn',
        n_jobs = -1,
        alpha = clamping_factor,
        n_neighbors = k
        )
        ls.fit(gex.values, partial_labels)
        res[:,i] = ls.transduction_
    full_labs = np.empty(res.shape[0])
    for i, row in enumerate(res):
        unique_elements, counts = np.unique(row, return_counts = True)
        recurring_class = unique_elements[np.argmax(counts)]
        full_labs[i] = recurring_class
    return pd.Series(full_labs, index = partial_labels.index)


def get_full_classif(gex, cells_by_fov, partial_labels, ks, clamping_factor = 0.2):
    p_res = []
    for fov, cells in cells_by_fov.items():
        p_gex = gex.loc[cells,:]
        p_labs = partial_labels.loc[cells]
        p_res.append(ensemble_ls(p_gex, p_labs, ks, clamping_factor))
    full_classif = pd.concat(p_res)
    return full_classif


def summarize_classification_result(classifications_dict, gsea_res_dict, confidence_scores_dict, classification_name):
    to_concat = []
    for smp in samples:
        df = pd.DataFrame({
            'sample_id' : smp,
            'cell_id' : classifications_dict[smp].index,
            'fov' : [x.split('_')[2] for x in classifications_dict[smp].index],
            classification_name : classifications_dict[smp].values
        })
        to_concat.append(df)
    res = pd.concat(to_concat, axis = 0)
    scores = pd.concat(confidence_scores_dict.values())
    info = pd.concat(gsea_res_dict.values(), axis = 1).T
    info.columns = [f'{st}_nes' for st in info.columns]
    info['confidence_score'] = scores
    info.reset_index(drop = False, inplace = True)
    info.rename({'Name' : 'cell_id'}, axis = 1, inplace = True)
    classification_res = pd.merge(res, info, how = 'inner', on = 'cell_id')
    return classification_res