# COMMOT
Screening cell-cell communication in spatial transcriptomics via collective optimal transport 

![pytest](https://github.com/zcang/COMMOT/actions/workflows/python-package.yml/badge.svg)

## Installation
Install by `pip install COMMOT` or

1. Install the dependencies. \
   `pip install -r requirements.txt`
2. Install COMMOT by cd to this directory and \
   `pip install .`
3. [Optional] Install [tradeSeq](https://github.com/statOmics/tradeSeq) in R to analyze the CCC differentially expressed genes. \
   Currently, tradeSeq version 1.0.1 with R version 3.6.3 has been tested to work. \
   In order for the R-python interface to work properly, rpy2==3.4.2 and anndata2ri==1.0.6 should be installed.

## Usage
**Basic usage**

_Import packages_
```
import commot as ct
import scanpy as sc
import pandas as pd
```
_Load a spatial dataset_ \
(e.g., a Visium dataset)
```
adata = sc.datasets.visium_sge(sample_id='V1_Breast_Cancer_Block_A_Section_1')
```
_Basic processing_
```
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
```
_Specify ligand-receptor pairs_
```
LR=np.array([['AMH', 'ACVR1'],['AMH', 'AMHR2'],['BMP10', 'ACVR1']],dtype=str)
df_ligrec = pd.DataFrame(data=LR)
```
(or use pairs from a ligand-receptor database `df_ligrec=ct.pp.ligand_receptor_database(database='CellChat', species='human')`.)

_Construct CCC networks_ \
Use collective optimal transport to construct CCC networks for the ligand-receptor pairs with a spatial distance constraint of 1000 (coupling between cells with distance greater than 1000 is prohibited). For example, the spot-by-spot matrix for the pair AMH (ligand) and ACVR1 (receptor)is stored in `adata.obsp['commot-user_database-AMH-ACVR1']`. The total sent or received signal for each pair is stored in `adata.obsm['commot-user_database-sum-sender']` and `adata.obsm['commot-user_database-sum-receiver']`.
```
ct.tl.spatial_communication(adata,
    pathway_name='user_database', df_ligrec=df_ligrec, dis_thr=1000)
```
**Documentation**

See the documentation in `docs/_build/html/index.html` for all the APIs to perform visualization and analyses such as visualizing spatial signaling direction and identifying CCC differentially expressed genes.

**Reference**
Cang, Zixuan, Yanxiang Zhao, Axel A. Almet, Adam Stabell, Raul Ramos, Maksim Plikus, Scott X. Atwood, and Qing Nie. "Screening cell-cell communication in spatial transcriptomics via collective optimal transport." bioRxiv (2022): 505185
