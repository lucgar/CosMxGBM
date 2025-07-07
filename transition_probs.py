import os
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import anndata as ad
import scanpy as sc
import squidpy as sq
import cellrank as cr


def calculate_transition_probs(adata):
    
    source_labs = [
        'clusteredGPM',
        'dispersedGPM',
        'clusteredMTC',
        'dispersedMTC',
        'clusteredNEU',
        'dispersedNEU',
        'clusteredPPR',
        'dispersedPPR',
    ]

    target_labs = [
        'GPM',
        'MTC',
        'NEU',
        'PPR'
    ]
    
    #copying adata
    tmp = ad.AnnData(X = adata.layers['raw_counts'], obs = adata.obs, var = adata.var)
    
    #filtering
    tmp = tmp[
        (tmp.obs.tumor_classification.isin(target_labs)) & 
        (~tmp.obs.struct_label.isna()) & 
        (tmp.obs.struct_label != 'other')
    ].copy()
    
    #prep
    sc.pp.normalize_total(tmp, target_sum = 1e5)
    print('finished norm')
    rsc.pp.pca(tmp, random_state = 37)
    print('finished pca')
    rsc.pp.neighbors(tmp, random_state = 37)
    print('finished neighs')
    rsc.tl.diffmap(tmp)
    print('finished diffmaps')
    rsc.pp.neighbors(tmp, use_rep = 'X_diffmap', key_added = 'dm', random_state = 37)
    print('finished diffmaps neighs')
    tmp.uns['iroot'] = int(np.where(tmp.obs.complexity == tmp.obs.complexity.max())[0][0])
    sc.tl.dpt(tmp, neighbors_key = 'dm')
    print('finished pseudotime')
    
    #calc pseudotime kernel
    pk = cr.kernels.PseudotimeKernel(tmp, time_key = 'dpt_pseudotime')
    pk.compute_transition_matrix()
    print(pk)
    
    #get transition matrix
    transition_matrix = pd.DataFrame(
        pk.transition_matrix.toarray(), 
        index = tmp.obs.struct_label.values, 
        columns = tmp.obs.tumor_classification.values
    )
    
    tmp_res = {}

    for sl in source_labs:
        tmp = transition_matrix.loc[sl]
        to_concat = []
        for tl in target_labs:
            pp = tmp.loc[:,tl].sum(axis = 1)#.mean()
            to_concat.append(pp)
        tt = pd.concat(to_concat, axis = 1)
        tt.columns = target_labs
        tmp_res[sl] = tt
    
    return tmp_res


# summarize transtions from result dict

source_labs = [
    'clusteredGPM',
    'dispersedGPM',
    'clusteredMTC',
    'dispersedMTC',
    'clusteredNEU',
    'dispersedNEU',
    'clusteredPPR',
    'dispersedPPR',
]

target_labs = [
    'GPM',
    'MTC',
    'NEU',
    'PPR'
]


rows = []
for sl in source_labs:
    to_concat = [results[smp][sl] for smp in results.keys()]
    tmp = pd.concat(to_concat, axis = 0)
    row = tmp.mean(axis = 0).values
    prob = row/row.sum()
    rows.append(row)
    
res = pd.DataFrame(rows, index = source_labs, columns = target_labs)