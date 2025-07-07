import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq


def calc_swed(adata, alpha = .3, random_state = 37):
    
    exp_graph = adata.obsp['connectivities']
    
    sq.gr.spatial_neighbors(
        adata,
        spatial_key = 'global_coords',
        coord_type = 'generic',
        n_neighs = 20
    )
    
    spat_graph = adata.obsp['spatial_connectivities']
    
    joint_graph = (1-alpha) * exp_graph + alpha * spat_graph
    
    sc.tl.louvain(adata, random_state = random_state, adjacency = joint_graph, use_weights = True, key_added = 'swed')