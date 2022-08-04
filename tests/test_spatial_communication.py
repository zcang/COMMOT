import pytest
import commot as ct

import anndata
import numpy as np
import scanpy as sc
import pandas as pd


def test_spatial_communication():

    X = np.array([[0.1, 0.2, 0.5, 1.0],
              [0.5, 0.2, 0.4, 1.0],
              [1.5, 0.2, 0.7, 0.0],
              [0.0, 0.0, 0.0, 0.0],
              [0.0, 1.0, 1.0, 0.0],
              [0.0, 1.0, 0.0, 0.0],
              [1.0, 0.1, 1.0, 0.1]], float)

    pts = np.array([[0.0, 0.0],
                    [0.0, 1.0],
                    [0.0, 2.0],
                    [1.0, 0.5],
                    [1.0, 0.0],
                    [1.2, 2.3],
                    [1.2, 0.7]], float)

    adata = anndata.AnnData(X=X, var=pd.DataFrame(index=['Lig1','Lig2','Rec1','Rec2']) )
    adata.obsm['spatial'] = pts
    ligrec =[['Lig1','Rec1','pathway1'],
            ['Lig1','Rec2','pathway1'],
            ['Lig2','Rec2','pathway2']]
    df_ligrec = pd.DataFrame(data=ligrec, columns=['ligand','receptor','pathway'])

    ct.tl.spatial_communication(adata, df_ligrec = df_ligrec, dis_thr=1.0, heteromeric=False, database_name='testdb', cot_eps_p=0.1, cot_nitermax=10000)

    X1 = [[0.00833303708452515, 2.32396343012789e-08, 0.0, 0.0, 0.07500079408819754, 0.0, 0.0], 
          [0.375026335084811, 0.035699164267244066, 2.2532037812087345e-10, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.3571357170372008, 0.699998961398071, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.01979814863634195, 0.0, 0.9576201052208858]]

    X2 = [[0.06666547680166099, 1.6269908894019704e-07, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.12511428752995418, 0.21414874159510383, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.5591574484245683, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.07030451439328901]]

    X3 = [[0.050023764100675105, 0.1499762034541506, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [2.461345888519154e-05, 0.19996577591605905, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.15883893013856648, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.9248906001487338, 0.0, 0.0, 0.0, 0.0, 0.0, 1.4275138992249157e-05], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0772631572277707]]

    X1 = np.array(X1)
    X2 = np.array(X2)
    X3 = np.array(X3)

    assert np.sum(np.abs(adata.obsp['commot-testdb-Lig1-Rec1'].toarray()-X1)) < 1e-5, "OT plan does not agree with benchmark"
    assert np.sum(np.abs(adata.obsp['commot-testdb-Lig1-Rec2'].toarray()-X2)) < 1e-5, "OT plan does not agree with benchmark"
    assert np.sum(np.abs(adata.obsp['commot-testdb-Lig2-Rec2'].toarray()-X3)) < 1e-5, "OT plan does not agree with benchmark"


