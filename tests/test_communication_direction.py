import pytest
import commot as ct

import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from scipy import sparse


def test_communication_direction():

    obsp_X = np.array([[0.86, 0.  , 0.91, 0.  , 0.  , 0.53, 0.  , 0.83, 0.  , 0.  ],
        [0.  , 0.  , 0.  , 0.97, 0.  , 0.  , 0.  , 0.82, 0.  , 0.  ],
        [0.  , 0.9 , 0.  , 0.99, 0.56, 0.  , 0.99, 0.  , 0.  , 0.  ],
        [0.55, 0.  , 0.74, 0.81, 0.  , 0.  , 0.  , 0.61, 0.75, 0.  ],
        [0.  , 0.  , 0.  , 0.  , 0.8 , 0.  , 0.86, 0.  , 0.91, 0.  ],
        [0.92, 0.94, 0.  , 0.67, 0.9 , 0.61, 0.54, 0.  , 0.98, 0.  ],
        [0.  , 0.56, 0.  , 0.  , 0.67, 0.  , 0.  , 0.  , 0.9 , 0.  ],
        [0.  , 0.66, 0.59, 0.72, 0.  , 0.86, 0.97, 0.  , 0.  , 0.81],
        [0.  , 0.  , 0.  , 0.  , 0.67, 0.  , 0.  , 0.  , 0.54, 0.96],
        [0.  , 0.85, 0.81, 0.  , 0.  , 0.  , 0.  , 0.98, 0.82, 0.  ]])

    pts = np.array([[0.98, 0.92],
        [0.67, 0.85],
        [0.31, 0.63],
        [0.56, 0.01],
        [0.51, 0.4 ],
        [0.42, 0.45],
        [0.29, 0.71],
        [0.85, 0.06],
        [0.88, 0.67],
        [0.73, 0.58]])
    
    vf_sender = np.array([[-2.87247023, -1.24330801],
        [ 0.05710527, -1.78908887],
        [ 3.27166403,  1.06292731],
        [ 1.50950633,  3.11335681],
        [ 0.48403577,  2.52400661],
        [ 4.16868451,  3.67908541],
        [ 1.9259737 , -0.90968418],
        [-3.20569482,  3.31294743],
        [-1.81820893, -1.18448988],
        [ 0.95583767, -3.32535327]])

    vf_receiver = np.array([[ 1.78472001,  1.49789001],
        [ 2.76449191,  2.76508308],
        [-3.01071112, -0.4879739 ],
        [ 0.53604126, -4.12531935],
        [ 3.14696619, -1.74831455],
        [-1.48143807,  1.34362988],
        [-1.53858178,  2.98703299],
        [ 0.175478  , -3.23524458],
        [ 4.21573862,  2.49750833],
        [-1.69823424,  0.49889926]])

    adata = anndata.AnnData(X = np.random.rand(10,2))
    adata.obsm['spatial'] = pts
    adata.obsp['commot-databaseX-total-total'] = sparse.csr_matrix(obsp_X)

    ct.tl.communication_direction(adata, database_name='databaseX', k=2)

    assert np.sum(np.abs(adata.obsm['commot_sender_vf-databaseX-total-total']-vf_sender)) < 1e-7, "Communication direction interpolation wrong."

    assert np.sum(np.abs(adata.obsm['commot_receiver_vf-databaseX-total-total']-vf_receiver)) < 1e-7, "Communication direction interpolation wrong."


