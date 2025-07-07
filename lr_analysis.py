import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import rapids_singlecell as rsc

res = rsc.gr.ligrec(
        idata,
        cluster_key = 'test_var',
        clusters = [
            ('GPMclustered', 'GPMclustered'),
            ('PPRclustered', 'PPRclustered'),
            ('dispersed', 'dispersed'),
        ],
        interactions = ligrec,
        threshold = 0,
        corr_method = 'fdr_bh',
        corr_axis = 'clusters',
        use_raw = False,
        n_perms = 1_000,
        copy = True
    )