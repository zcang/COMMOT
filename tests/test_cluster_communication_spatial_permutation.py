import pytest
import commot as ct

import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from scipy import sparse


def test_cluster_communication():

    X = np.array([[0.24, 0.9 , 0.44, 0.99],
        [0.56, 0.2 , 0.99, 0.28],
        [0.34, 0.08, 0.55, 0.44],
        [0.74, 0.81, 0.49, 0.4 ],
        [0.07, 0.61, 0.75, 0.26],
        [0.37, 0.44, 0.12, 0.01],
        [0.8 , 0.19, 0.86, 0.47],
        [0.91, 0.27, 0.92, 0.94],
        [0.38, 0.67, 0.9 , 0.61],
        [0.54, 0.2 , 0.98, 0.24]])
    
    pts = np.array([[0.86, 0.26],
        [0.91, 0.  ],
        [0.31, 0.53],
        [0.37, 0.83],
        [0.31, 0.16],
        [0.4 , 0.19],
        [0.46, 0.97],
        [0.39, 0.2 ],
        [0.17, 0.82],
        [0.2 , 0.18]])
        
    comm_mat_1 = np.array([[0.12284352, 0.01066791, 0.03998729],
        [0.01884018, 0.12542461, 0.04087296],
        [0.01200533, 0.0525465 , 0.21689502]])
    
    pval_mat_1 = np.array(
        [[0.66, 0.75, 0.23],
        [0.56, 0.69, 0.49],
        [0.7, 0.41, 0.06]])
    
    comm_mat_2 = np.array([[0.17610059, 0.01235253, 0.05327244],
        [0.03120265, 0.16598969, 0.06185901],
        [0.02305418, 0.06354936, 0.30368174]])

    pval_mat_2 = np.array(
        [[0.64, 0.83, 0.33],
        [0.5,  0.77, 0.47],
        [0.65, 0.49, 0.08]])

    comm_mat_3 = np.array([[0.12284352, 0.01066791, 0.03998729],
        [0.01884018, 0.12542461, 0.04087296],
        [0.01200533, 0.0525465 , 0.21689502]])

    pval_mat_3 = np.array(
        [[0.36, 0.93, 0.21],
        [0.85, 0.19, 0.79],
        [0.66, 0.3,  0.33]])

    comm_mat_4 = np.array([[0.17610059, 0.01235253, 0.05327244],
        [0.03120265, 0.16598969, 0.06185901],
        [0.02305418, 0.06354936, 0.30368174]])

    pval_mat_4 = np.array(
        [[0.4,  0.94, 0.27],
        [0.89, 0.2,  0.73],
        [0.69, 0.28, 0.43]])

    celltype_id = np.array( [1,1,1,1,2,2,2,3,3,3], str )

    adata = anndata.AnnData(X=X, var=pd.DataFrame(index=['Lig1','Lig2','Rec1','Rec2']) )
    adata.obsm['spatial'] = pts
    ligrec =[['Lig1','Rec1','pathway1'],
            ['Lig1','Rec2','pathway1'],
            ['Lig2','Rec2','pathway2']]
    df_ligrec = pd.DataFrame(data=ligrec, columns=['ligand','receptor','pathway'])
    adata.obs['celltype'] = pd.Series(list(celltype_id), dtype="category").values

    ct.tl.cluster_communication_spatial_permutation(adata, df_ligrec=df_ligrec, dis_thr=0.5, 
        database_name='databaseX', perm_type='all_cell', clustering='celltype', random_seed=1)

    assert np.sum(np.abs(adata.uns['commot_cluster_spatial_permutation-celltype-databaseX-Lig1-Rec1']['communication_matrix'].values-comm_mat_1)) < 1e-7, "Cluster communication score wrong (all_cell permutaton)."
    assert np.sum(np.abs(adata.uns['commot_cluster_spatial_permutation-celltype-databaseX-pathway1']['communication_matrix'].values-comm_mat_2)) < 1e-7, "Cluster communication score wrong (all_cell permutaton)."
    assert np.sum(np.abs(adata.uns['commot_cluster_spatial_permutation-celltype-databaseX-Lig1-Rec1']['communication_pvalue'].values-pval_mat_1)) < 1e-7, "Cluster communication p-value wrong (all_cell permutaton)."
    assert np.sum(np.abs(adata.uns['commot_cluster_spatial_permutation-celltype-databaseX-pathway1']['communication_pvalue'].values-pval_mat_2)) < 1e-7, "Cluster communication p-value wrong (all_cell permutaton)."

    ct.tl.cluster_communication_spatial_permutation(adata, df_ligrec=df_ligrec, dis_thr=0.5, 
        database_name='databaseX', perm_type='within_cluster', clustering='celltype', random_seed=1)

    assert np.sum(np.abs(adata.uns['commot_cluster_spatial_permutation-celltype-databaseX-Lig1-Rec1']['communication_matrix'].values-comm_mat_3)) < 1e-7, "Cluster communication score wrong (within_cluster permutaton)."
    assert np.sum(np.abs(adata.uns['commot_cluster_spatial_permutation-celltype-databaseX-pathway1']['communication_matrix'].values-comm_mat_4)) < 1e-7, "Cluster communication score wrong (within_cluster permutaton)."
    assert np.sum(np.abs(adata.uns['commot_cluster_spatial_permutation-celltype-databaseX-Lig1-Rec1']['communication_pvalue'].values-pval_mat_3)) < 1e-7, "Cluster communication p-value wrong (within_cluster permutaton)."
    assert np.sum(np.abs(adata.uns['commot_cluster_spatial_permutation-celltype-databaseX-pathway1']['communication_pvalue'].values-pval_mat_4)) < 1e-7, "Cluster communication p-value wrong (within_cluster permutaton)."
