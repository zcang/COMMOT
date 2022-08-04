import pytest
import commot as ct

import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from scipy import sparse


def test_cluster_communication():

    obsp_X = np.array([[0.  , 0.  , 0.  , 0.  , 0.  , 0.53, 0.37, 0.83, 0.31, 0.16],
        [0.  , 0.  , 0.  , 0.  , 0.  , 0.2 , 0.17, 0.82, 0.2 , 0.18],
        [0.  , 0.  , 0.  , 0.  , 0.  , 0.2 , 0.99, 0.28, 0.34, 0.08],
        [0.  , 0.  , 0.  , 0.  , 0.  , 0.4 , 0.07, 0.61, 0.75, 0.26],
        [0.  , 0.  , 0.  , 0.  , 0.  , 0.19, 0.86, 0.47, 0.91, 0.27],
        [0.  , 0.  , 0.  , 0.  , 0.  , 0.61, 0.54, 0.2 , 0.98, 0.24],
        [0.  , 0.  , 0.  , 0.  , 0.  , 0.11, 0.34, 0.38, 0.9 , 0.41],
        [0.28, 0.66, 0.59, 0.72, 0.22, 0.  , 0.  , 0.  , 0.  , 0.  ],
        [0.44, 0.16, 0.23, 0.37, 0.67, 0.  , 0.  , 0.  , 0.  , 0.  ],
        [0.14, 0.85, 0.81, 0.02, 0.26, 0.  , 0.  , 0.  , 0.  , 0.  ]])

    comm_mat = np.array([[0.        , 0.24416667, 0.40166667],
        [0.        , 0.29444444, 0.52888889],
        [0.43916667, 0.12777778, 0.        ]])

    pval_mat = np.array([[1.  , 0.42, 0.04],
        [1.  , 0.05, 0.02],
        [0.  , 0.91, 1.  ]])

    celltype_id = np.array( [1,1,1,1,2,2,2,3,3,3], str )
    adata = anndata.AnnData(X=np.random.rand(10,2))
    adata.obsp['commot-databaseX-total-total'] = sparse.csr_matrix(obsp_X)
    adata.obs['celltype'] = pd.Series(list(celltype_id), dtype="category").values

    ct.tl.cluster_communication(adata, clustering='celltype', random_seed=1, database_name='databaseX', n_permutations = 100)

    assert np.sum(np.abs(adata.uns['commot_cluster-celltype-databaseX-total-total']['communication_matrix'].values-comm_mat)) < 1e-7, "Cluster communication score wrong."

    assert np.sum(np.abs(adata.uns['commot_cluster-celltype-databaseX-total-total']['communication_pvalue'].values-pval_mat)) < 1e-7, "Cluster commuinication p-value wrong."
