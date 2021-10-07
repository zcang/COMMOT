# COMMOT
Cell-cell communication network construction for spatial transcriptomics using collective optimal transport

## Installation
1. Install the dependencies. \
   `pip install -r requirements.txt`
2. Install COMMOT by cd to this directory and \
   `pip install .`
3. Install tradeSeq in R to analyze the CCC differentially expressed genes. \
   Currently, tradeSeq version 1.0.1 with R version 3.6.3 has been tested to work. \
   In order for the R-python interface to work properly, rpy2==3.4.2 and anndata2ri==1.0.6 should be installed.

## Usage
**Basic usage**

Import packages
```
import commot
import scanpy as sc
```
Load a spatial dataset (e.g., a Visium dataset)
```
adata = sc.datasets.visium_sge(sample_id='V1_Breast_Cancer_Block_A_Section_1')
```
