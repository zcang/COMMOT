# COMMOT
Screening cell-cell communication in spatial transcriptomics via collective optimal transport 

![pytest](https://github.com/zcang/COMMOT/actions/workflows/python-package.yml/badge.svg)

## Installation
Install from PyPI by `pip install commot` or install from source by cd to this directory and `pip install .`

[Optional] Install [tradeSeq](https://github.com/statOmics/tradeSeq) in R to analyze the CCC differentially expressed genes. \
Currently, tradeSeq version 1.0.1 with R version 3.6.3 has been tested to work. \
In order for the R-python interface to work properly, rpy2==3.4.2 and anndata2ri==1.0.6 should be installed.
To use this feature, install from PyPI by `pip install commot[tradeSeq]` or from source by `pip install .[tradeSeq]`.

## Usage
**Basic usage**

_Import packages_
```
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
```
_Load a spatial dataset_ \
(e.g., a Visium dataset)
```
adata = sc.datasets.visium_sge(sample_id='V1_Mouse_Brain_Sagittal_Posterior')
adata.var_names_make_unique()
```
_Basic processing_
```
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
```
_Specify ligand-receptor pairs_
```
LR=np.array([['Tgfb1', 'Tgfbr1_Tgfbr2', 'Tgfb_pathway'],['Tgfb2', 'Tgfbr1_Tgfbr2', 'Tgfb_pathway'],['Tgfb3', 'Tgfbr1_Tgfbr2', 'Tgfb_pathway']],dtype=str)
df_ligrec = pd.DataFrame(data=LR)
```
(or use pairs from a ligand-receptor database `df_ligrec=ct.pp.ligand_receptor_database(database='CellChat', species='mouse')`.)

_Construct CCC networks_ \
Use collective optimal transport to construct CCC networks for the ligand-receptor pairs with a spatial distance constraint of 200 (coupling between cells with distance greater than 200 is prohibited). For example, the spot-by-spot matrix for the pair Tgfb1 (ligand) and Tgfbr1_Tgfbr2 (receptor)is stored in `adata.obsp['commot-user_database-Tgfb1-Tgfbr1_Tgfbr2']`. The total sent or received signal for each pair is stored in `adata.obsm['commot-user_database-sum-sender']` and `adata.obsm['commot-user_database-sum-receiver']`.
```
ct.tl.spatial_communication(adata,
    database_name='user_database', df_ligrec=df_ligrec, dis_thr=200, heteromeric=True)
```
**Documentation**

See the documentation at [https://commot.readthedocs.io/en/latest/index.html](https://commot.readthedocs.io/en/latest/index.html) for all the APIs to perform visualization and analyses such as visualizing spatial signaling direction and identifying CCC differentially expressed genes.

**Reference**

Cang, Zixuan, Yanxiang Zhao, Axel A. Almet, Adam Stabell, Raul Ramos, Maksim Plikus, Scott X. Atwood, and Qing Nie. "Screening cell-cell communication in spatial transcriptomics via collective optimal transport." bioRxiv (2022): 505185
