# Installation

## Basic installation
Python 3.7 (or later) is required. Consider using virtual environments with [Miniconda] or [Anaconda].

To install, run
```
pip install commot
```
Alternatively, you can install the latest version from [GitHub]
```
git clone https://github.com/zcang/COMMOT.git
cd COMMOT
pip install .
```
## Optional dependency
One of the downstream analysis function `commot.tl.communication_deg_detection()` relies on tradeSeq [[1]](#1) to find signaling DE genes.
To use this function, tradeSeq and some extra python packages need to be installed.

To install commot with extra packages, run
```
pip install commot[tradeSeq]
```

Currently, the R-python interface was tested to only work with R-3.6.3 and tradeSeq-1.0.1. Consider the following approach to install R-3.6.3 together with other R versions, for example, on Ubuntu 20.04, run
```
export R_VERSION=3.6.3
curl -O https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-${R_VERSION}_1_amd64.deb
sudo gdebi r-${R_VERSION}_1_amd64.deb
sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript
```
To switch between R versions, export `R_VERSION` to the desired version, remove the symbolic links `/usr/local/bin/R` and `/usr/local/bin/Rscript`, and run the last two lines again. Once R-3.6.3 is installed, install tradeSeq-1.0.1 by running
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tradeSeq")
```
We are working on an easier solution of this.

## Installation on Windows
**rvlib**\
`rvlib` may need to be installed prior. Install it by
```
conda install -c conda-forge rvlib
```

**pygeos**\
For Python 3.8+, `pygeos` may need to be uninstalled after installing `commot` due to incompatibilities with `Shapely`

**rasterio**\
For Python 3.7, `rasterio` may require installation using the appropriate whl files. See this [discussion](https://github.com/rasterio/rasterio/issues/1963#issuecomment-672262445) for installation instructions.

## References
<a id="1">[1]</a> 
Van den Berge, Koen, et al. "Trajectory-based differential expression analysis for single-cell sequencing data." Nature communications 11.1 (2020): 1-13.

[miniconda]: http://conda.pydata.org/miniconda.html
[anaconda]: https://www.anaconda.com/products/distribution
[github]: https://github.com/zcang/commot