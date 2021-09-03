from typing import Optional, Union
import ot
import sys
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
from matplotlib import cm
import matplotlib.pyplot as plt
import plotly
from scipy import sparse
from scipy.spatial import distance_matrix
from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import normalize
from sklearn.neighbors import kneighbors_graph
from sklearn.decomposition import PCA
import karateclub

from .._utils import partial_corr
from .._utils import semipartial_corr
from .._utils import treebased_score
from .._utils import treebased_score_multifeature
from .._utils import d_graph_local_jaccard
from .._utils import d_graph_local_jaccard_weighted
from .._utils import d_graph_global_structure
from .._utils import leiden_clustering
from .._utils import moranI_vector_global
from .._utils import preprocess_vector_field
from .._utils import binarize_sparse_matrix

def communication_deg_detection(
    adata: anndata.AnnData,
    n_var_genes: int = 3000,
    pathway_name: str = None,
    summary: str = 'receiver',
    lr_pair: tuple = ('total','total'),
    nknots: int = 6,
    n_deg_genes: int = None,
    n_points: int = 50,
    deg_pvalue_cutoff: float = 0.05,
):
    """
    Identify signaling dependent genes

    Paramters
    ---------
    adata
        The data matrix of shape ``n_obs`` × ``n_var``.
        Rows correspond to cells or positions and columns to genes.
        The count data should be available through adata.layers['count'].
        The signaling data should be available in 
        adata.obsm['commot-$pathway_name-sum']['$summary-$ligand-$receptor']
    n_var_genes
        The number of most variable genes to test.
    n_deg_genes
        The number of top deg genes to evaluate yhat.
    pathway_name
        Name of the signaling pathway.
    summary
        'sender' or 'receiver'
    lr_pair
        A tuple of the ligand-receptor pair
    nknots
        Number of knots in spline when constructing GAM
    n_points
        Number of points on which to evaluate the fitted GAM 
        for downstream clustering and visualization.
    deg_pvalue_cutoff
        The p-value cutoff of genes for obtaining the fitted gene expression patterns.

    Returns
    -------
    df_deg: pd.DataFrame
        A data frame of deg analysis results, including Wald statistics, degree of freedom, and p-value.
    df_yhat: pd.DataFrame
        A data frame of smoothed gene expression values.
    
    """
    # setup R environment
    # !!! anndata2ri works only with 3.6.3 in the tested machine
    import rpy2
    import anndata2ri
    import rpy2.robjects as ro
    from rpy2.robjects.conversion import localconverter
    import rpy2.rinterface_lib.callbacks
    import logging
    rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

    ro.r('library(tradeSeq)')
    ro.r('library(clusterExperiment)')
    anndata2ri.activate()
    ro.numpy2ri.activate()
    ro.pandas2ri.activate()

    # prepare input adata for R
    adata_deg = anndata.AnnData(
        X = adata.layers['counts'],
        var = pd.DataFrame(index=list(adata.var_names)),
        obs = pd.DataFrame(index=list(adata.obs_names)))
    adata_deg_var = adata_deg.copy()
    sc.pp.log1p(adata_deg_var)
    sc.pp.highly_variable_genes(adata_deg_var, n_top_genes=n_var_genes)
    adata_deg = adata_deg[:, adata_deg_var.var.highly_variable]
    del adata_deg_var

    comm_sum = adata.obsm['commot-%s-sum' % pathway_name]['%s-%s-%s' % (summary, lr_pair[0], lr_pair[1])].values.reshape(-1,1)
    cell_weight = np.ones_like(comm_sum).reshape(-1,1)

    # send adata to R
    adata_r = anndata2ri.py2rpy(adata_deg)
    ro.r.assign("adata", adata_r)
    ro.r("X <- as.matrix( assay( adata, 'X') )")
    ro.r.assign("pseudoTime", comm_sum)
    ro.r.assign("cellWeight", cell_weight)

    # perform analysis (tradeSeq-1.0.1 in R-3.6.3)
    string_fitGAM = 'sce <- fitGAM(counts=X, pseudotime=pseudoTime[,1], cellWeights=cellWeight[,1], nknots=%d, verbose=TRUE)' % nknots
    ro.r(string_fitGAM)
    ro.r('assoRes <- data.frame( associationTest(sce, global=FALSE, lineage=TRUE) )')
    ro.r('assoRes[is.nan(assoRes[,"waldStat_1"]),"waldStat_1"] <- 0.0')
    ro.r('assoRes[is.nan(assoRes[,"df_1"]),"df_1"] <- 0.0')
    ro.r('assoRes[is.nan(assoRes[,"pvalue_1"]),"pvalue_1"] <- 1.0')
    with localconverter(ro.pandas2ri.converter):
        df_assoRes = ro.r['assoRes']
    ro.r('assoRes = assoRes[assoRes[,"pvalue_1"] <= %f,]' % deg_pvalue_cutoff)
    ro.r('oAsso <- order(assoRes[,"waldStat_1"], decreasing=TRUE)')
    if n_deg_genes is None:
        n_deg_genes = df_assoRes.shape[0]
    string_cluster = 'clusPat <- clusterExpressionPatterns(sce, nPoints = %d,' % n_points\
        + 'verbose=TRUE, genes = rownames(assoRes)[oAsso][1:min(%d,length(oAsso))],' % n_deg_genes \
        + ' k0s=4:5, alphas=c(0.1))'
    ro.r(string_cluster)
    ro.r('yhatScaled <- data.frame(clusPat$yhatScaled)')
    with localconverter(ro.pandas2ri.converter):
        yhat_scaled = ro.r['yhatScaled']

    df_deg = df_assoRes.rename(columns={'waldStat_1':'waldStat', 'df_1':'df', 'pvalue_1':'pvalue'})
    idx = np.argsort(-df_deg['waldStat'].values)
    df_deg = df_deg.iloc[idx]
    df_yhat = yhat_scaled

    anndata2ri.deactivate()
    ro.numpy2ri.deactivate()
    ro.pandas2ri.deactivate()

    return df_deg, df_yhat
    
def communication_deg_clustering(
    df_deg: pd.DataFrame,
    df_yhat: pd.DataFrame,
    deg_clustering_npc: int = 10,
    deg_clustering_knn: int = 5,
    deg_clustering_res: float = 1.0,
    n_deg_genes: int = 200,
    p_value_cutoff: float = 0.05
):
    """
    Cluster the communcation DE genes based on their fitted expression pattern.

    Parameters
    ----------
    df_deg
        The deg analysis summary data frame obtained by running ``communication_deg_detection``.
        Each row corresponds to one tested genes and columns include "waldStat" (Wald statistics), "df" (degrees of freedom), and "pvalue" (p-value of the Wald statistics).
    df_yhat
        The fitted (smoothed) gene expression pattern obtained by running ``tl.communication_deg_detection``.
    deg_clustering_npc
        Number of PCs when performing PCA to cluster gene expression patterns
    deg_clustering_knn
        Number of neighbors when constructing the knn graph for leiden clustering.
    deg_clustering_res
        The resolution parameter for leiden clustering.
    n_deg_genes
        Number of top deg genes to cluster.
    p_value_cutoff
        The p-value cutoff for genes to be included in clustering analysis.

    Returns
    -------
    df_deg_clus: pd.DataFrame
        A data frame of clustered genes.
    df_yhat_clus: pd.DataFrame
        The fitted gene expression patterns of the clustered genes

    """
    df_deg = df_deg[df_deg['pvalue'] <= p_value_cutoff]
    n_deg_genes = min(n_deg_genes, df_deg.shape[0])
    idx = np.argsort(-df_deg['waldStat'])
    df_deg = df_deg.iloc[idx[:n_deg_genes]]
    yhat_scaled = df_yhat.loc[df_deg.index]
    x_pca = PCA(n_components=deg_clustering_npc, svd_solver='full').fit_transform(yhat_scaled.values)
    d_pca = distance_matrix(x_pca, x_pca)
    cluster_labels = leiden_clustering(d_pca, k=deg_clustering_knn, resolution=deg_clustering_res)

    data_tmp = np.concatenate((df_deg.values, cluster_labels.reshape(-1,1)),axis=1)
    df_metadata = pd.DataFrame(data=data_tmp, index=df_deg.index,
        columns=['waldStat','df','pvalue','cluster'] )
    return df_metadata, yhat_scaled

def communication_impact(
    adata: anndata.AnnData,
    pathway_name: str = None,
    normalize: bool=False,
    method: str = None,
    corr_method: str = "spearman",
    tree_method: str = "rf",
    tree_ntrees: int = 100,
    tree_repeat: int = 100,
    tree_max_depth: int = 5,
    tree_max_features: str = 'sqrt',
    tree_learning_rate: float = 0.1,
    tree_subsample: float = 1.0,
    tree_combined: bool = False,
    ds_genes: list = None,
    bg_genes: Union[list, int] = 100
):
    """
    Analyze impact of communication

    Parameters
    ----------
    adata
        The data matrix of shape ``n_obs`` × ``n_var``.
        Rows correspond to cells or positions and columns to genes.
        The full normalized dataset should be available in ``adata.raw``.
    pathway_name
        Name of the signaling pathway.
    normalize
        Whether to perform normalization before determining variable genes.
    method
        'partial_corr': partial correlation.
        'semipartial_corr': semipartial correlation.
        'treebased_score': machine learning based score (ensemble of trees).
    corr_method
        The correlation coefficient to use when method is 'partial_corr' or 'semipartial_corr'.
        'spearman': Spearman's r. 'pearson': Pearson's r.
    tree_method
        The ensemble of trees method to use when method is 'treebased_score'.
        'gbt': gradient boosted trees. 'rf': random forest.
    tree_ntrees
        Number of trees when using 'treebased_score'.
    tree_repeat
        Number of times to repeat to account for randomness when using 'treebased_score'.
    tree_mas_depth
        Max depth of trees when using 'treebased_score'.
    tree_max_features
        Max features for trees when using 'treebased_score'.
    tree_learning_rate
        Learning rate when using 'treebased_score'.
    tree_subsample
        Subsample (between 0 and 1) when using 'treebased_score'.
    tree_combined
        If True, use a single model for each target gene with all features.
    ds_genes
        A list of genes for analyzing the correlation with cell-cell communication. 
    bg_genes
        If an integer, the top number of variable genes are used.
        Alternatively, a list of genes.

    Returns
    -------
    df_impact: pd.DataFrame
        A data frame describing the correlation 
        between the ds_genes and cell-cell communication.
    """

    # Get a list of background genes using most 
    # variable genes if only given a number.
    adata_bg = adata.raw.to_adata()
    adata_all = adata.raw.to_adata()
    if normalize:
        sc.pp.normalize_total(adata_bg, inplace=True)
        sc.pp.log1p(adata_bg)
    if np.isscalar(bg_genes):
        ng_bg = int(bg_genes)
        sc.pp.highly_variable_genes(adata_bg, n_top_genes=ng_bg)
        adata_bg = adata_bg[:,adata_bg.var.highly_variable]
    else:
        adata_bg = adata_bg[:,bg_genes]
    # Prepare downstream or upstream genes
    ncell = adata.shape[0]
    col_names = []
    Ds_exps = []
    Ds_exp_total = np.zeros([ncell], float)
    for i in range(len(ds_genes)):
        Ds_exp = np.array(adata_all[:,ds_genes[i]].X.toarray()).reshape(-1)
        Ds_exps.append(Ds_exp)
        col_names.append(ds_genes[i])
        Ds_exp_total += Ds_exp
    Ds_exps.append(Ds_exp_total); col_names.append('average')
    # Impact analysis
    bg_genes = list( adata_bg.var_names )
    nrows = adata.obsm['commot-'+pathway_name+'-sum'].shape[1]
    ncols = len(ds_genes) + 1
    impact_mat = np.empty([nrows, ncols], float)
    df_ligrec = adata.uns['commot-'+pathway_name+'-info']['df_ligrec']
    row_names_sender = []; row_names_receiver = []
    exclude_lr_genes_list = []
    for i in range(df_ligrec.shape[0]):
        lig = df_ligrec.iloc[i][0]
        rec = df_ligrec.iloc[i][1]
        row_names_sender.append('sender-%s-%s' % (lig, rec))
        row_names_receiver.append('receiver-%s-%s' % (lig, rec))
        exclude_lr_genes_list.append(set((lig,rec)))
    row_names_sender.append('sender-total-total')
    row_names_receiver.append('receiver-total-total')
    exclude_lr_genes_list.append(set(tuple(list(df_ligrec.values.flatten()))))
    row_names = row_names_sender + row_names_receiver
    exclude_lr_genes_list = exclude_lr_genes_list + exclude_lr_genes_list


    for j in range(ncols):
        if j == ncols-1:
            exclude_ds_genes = set(ds_genes)
        else:
            exclude_ds_genes = set([ds_genes[j]])
        if method == 'treebased_score' and tree_combined:
            exclude_lr_genes = set(tuple(list(df_ligrec.values.flatten())))
            exclude_genes = list(exclude_lr_genes.union(exclude_ds_genes))
            use_genes = list( set(bg_genes) - set(exclude_genes) )
            bg_mat = np.array( adata_bg[:,use_genes].X.toarray() )
            r = treebased_score_multifeature(adata.obsm['commot-'+pathway_name+'-sum'][row_names].values, Ds_exps[j], bg_mat,
                n_trees = tree_ntrees, n_repeat = tree_repeat,
                max_depth = tree_max_depth, max_features = tree_max_features,
                learning_rate = tree_learning_rate, subsample = tree_subsample)
            impact_mat[:,j] = r[:]
        else:
            for i in range(nrows):
                row_name = row_names[i]
                exclude_lr_genes = exclude_lr_genes_list[i]

                exclude_genes = list(exclude_lr_genes.union(exclude_ds_genes))
                use_genes = list( set(bg_genes) - set(exclude_genes) )
                bg_mat = np.array( adata_bg[:,use_genes].X.toarray() )
                if method == "partial_corr":
                    r,p = partial_corr(adata.obsm['commot-'+pathway_name+'-sum'][row_name].values.reshape(-1,1), Ds_exps[j].reshape(-1,1), bg_mat, method=corr_method)
                elif method == "semipartial_corr":
                    r,p = semipartial_corr(adata.obsm['commot-'+pathway_name+'-sum'][row_name].values.reshape(-1,1), Ds_exps[j].reshape(-1,1), ycov=bg_mat, method=corr_method)
                elif method == "treebased_score":
                    r = treebased_score(adata.obsm['commot-'+pathway_name+'-sum'][row_name].values.reshape(-1,1), Ds_exps[j], bg_mat,
                        n_trees = tree_ntrees, n_repeat = tree_repeat,
                        max_depth = tree_max_depth, max_features = tree_max_features,
                        learning_rate = tree_learning_rate, subsample = tree_subsample)
                impact_mat[i,j] = r
    df_impact = pd.DataFrame(data=impact_mat, index = row_names, columns = col_names)
    return df_impact



def group_cluster_communication(
    adata: anndata.AnnData,
    pathway_name: Union[str, list] = None,
    clustering: str = None,
    keys = None,
    p_value_cutoff: float = 0.05,
    quantile_cutoff: float = 0.99,
    dissimilarity_method: str = None,
    leiden_k: int = 5,
    leiden_resolution: float = 1.0,
    leiden_random_seed: int = 1,
    leiden_n_iterations: int = -1,
    d_global_structure_weights: tuple = (0.45,0.45,0.1)
):
    """
    Idenfitify groups of cluster-cluster communication with similar
    pattern.

    Parameters
    ----------
    adata
        The data matrix with the cluster-cluster communication 
        info stored in ``adata.uns``.
    pathway_name
        A string for one signaling pathway or a list of pathways.
    clustering
        Name of clustering with the labels stored in ``.obs[clustering]``.
    keys
        A list of keys for the analyzed communication connections as tuples 
        (pathway_name, ligand, receptor). 
        If given, pathway_name will be ignored.
    quantile_cutoff
        The quantile cutoff for including an edge. Set to 1 to disable this criterion.
        The quantile_cutoff and p_value_cutoff works in the "or" logic to avoid missing
        significant signaling connections.
    p_value_cutoff
        The cutoff of p-value to plot an edge.
    dissimilarity_method
        The dissimilarity measurement between graphs to use. 
        'jaccard' for Jaccard distance.
        'jaccard_weighted' for weighted Jaccard distance.
        'global_structure' for a metric focusing on global structure [Schieber2017]_.
    leiden_k
        Number of neighbors for the knn-graph to be fed to leiden clustering algorithm.
    leiden_resolution
        The resolution parameter for the leiden clustering algorithm.
    leiden_random_seed
        The random seed for the leiden clustering algorithm.
    leiden_n_iterations
        The maximum number of iterations for the leiden algorithm.
        The algorithm will run until convergence if set to -1.
    d_global_structure_weights
        The weights for the three terms in the global structural dissimilarity.
        See [Schieber2017]_ for more information.
    
    Returns
    -------
    keys : list
        The list of keys for the analyzed communication connections as tuples 
        (pathway_name, ligand, receptor).
    communication_clusterid : np.ndarray
        The group id of the cluster-cluster communications.
    D : np.ndarray
        The dissimilarity matrix for the cluster-cluster communications.

    References
    ----------
    .. [Schieber2017] Schieber, T. A., Carpi, L., Díaz-Guilera, A., Pardalos, 
        P. M., Masoller, C., & Ravetti, M. G. (2017). Quantification 
        of network structural dissimilarities. Nature communications, 8(1), 1-10.

    """
    # Get a list of keys if only pathway_name is specified.
    if keys is None:
        keys = []
        if isinstance(pathway_name, str):
            pathways = [pathway_name]
        elif isinstance(pathway_name, list):
            pathways = pathway_name
        for pathway in pathways:
            df_ligrec = adata.uns['commot-%s-info' % pathway]['df_ligrec']
            for i in range(df_ligrec.shape[0]):
                key = (pathway, df_ligrec.iloc[i][0], df_ligrec.iloc[i][1])
                keys.append(key)
            keys.append( (pathway,'total','total') )
    # Get a list of filtered communication matrices corresponding to the keys.
    As = []
    for key in keys:
        tmp_name = 'commot_cluster-%s-%s-%s-%s' % (key[0],clustering,key[1],key[2])
        X_tmp = adata.uns[tmp_name]['communication_matrix'].values.copy()
        pvalue_tmp = adata.uns[tmp_name]['communication_pvalue'].values.copy()
        if not quantile_cutoff is None:
            cutoff = np.quantile(X_tmp.reshape(-1), quantile_cutoff)
        else:
            cutoff = np.inf
        tmp_mask = ( X_tmp < cutoff ) * ( pvalue_tmp > p_value_cutoff )
        X_tmp[tmp_mask] = 0
        As.append(X_tmp)
    # Get a distance/dissimilarity matrix between the communication matrices.
    D = np.zeros([len(keys), len(keys)], float)
    for i in range(len(keys)-1):
        for j in range(i+1,len(keys)):
            if dissimilarity_method == 'jaccard':
                d = d_graph_local_jaccard(As[i], As[j])
            elif dissimilarity_method == 'jaccard_weighted':
                d = d_graph_local_jaccard_weighted(As[i], As[j])
            elif dissimilarity_method == 'global_structure':
                w1,w2,w3 = d_global_structure_weights
                d = d_graph_global_structure(As[i], As[j], w1=w1, w2=w2, w3=w3)
            D[i,j] = d; D[j,i] = d
    # Perform clustering
    leiden_k = min(leiden_k, len(keys)-1)
    communication_clusterid = leiden_clustering(D,
        k = leiden_k, resolution = leiden_resolution,
        random_seed = leiden_random_seed, 
        n_iterations = leiden_n_iterations)
    
    return keys, communication_clusterid, D

def group_cell_communication(
    adata: anndata.AnnData,
    pathway_name: Union[str, list] = None,
    keys = None,
    bin_method: str = 'gaussian_mixture',
    bin_append_zeros: str = 'full',
    bin_random_state: int = 1,
    bin_cutoff: float = 0,
    knn: int = 2,
    dissimilarity_method: str = 'graphwave',
    kw_graphwave: dict = {},
    leiden_k: int = 5,
    leiden_resolution: float = 1.0,
    leiden_random_seed: int = 1,
    leiden_n_iterations: int = -1
):
    """
    Idenfitify groups of cell-cell communication with similar
    pattern.

    Parameters
    ----------
    adata
        The data matrix with the cell-cell communication 
        info stored in ``adata.obsp``.
    pathway_name
        A string for one signaling pathway or a list of pathways.
    keys
        A list of keys for the analyzed communication connections as tuples 
        (pathway_name, ligand, receptor). 
        If given, pathway_name will be ignored.
    bin_method
        Method for binarize communication connections. Choices: 'gaussian_mixture', 'kmeans'.
    bin_append_zeros
        Number of zeros to append to the non-zero entries when running the binarization.
        'full' use the full flattened cell-by-cell communication matrix.
        'match' append the same number of zeros to the vector of non-zero entries.
    bin_random_state
        The random seed for binarization methods.
    bin_cutoff
        Force all connections with a weight smaller than or equal to bin_cutoff to be zero
        , regardless of binarization result.
    knn
        Number of neighbors for building the spatial knn graph.
    dissimilarity_method
        The method for quantifying dissimilarity.
        'graphwave', node structural embedding by GraphWave [Donnat2018]_ implemented 
        in `Karate Club
        <https://github.com/benedekrozemberczki/karateclub/>`_.
    kw_graphwave
        Keywords for GraphWave. Defaults: {'sample_number':200, 'step_size':0.1, 'heat_coefficient': 1.0
        'approximation':100, 'mechanism':'approximation', 'switch':1000, 'seed':42} See details at `Karate Club
        <https://github.com/benedekrozemberczki/karateclub/>`_.
    leiden_k
        Number of neighbors for the knn-graph to be fed to leiden clustering algorithm.
    leiden_resolution
        The resolution parameter for the leiden clustering algorithm.
    leiden_random_seed
        The random seed for the leiden clustering algorithm.
    leiden_n_iterations
        The maximum number of iterations for the leiden algorithm.
        The algorithm will run until convergence if set to -1.
    
    Returns
    -------
    keys : list
        The list of keys for the analyzed communication connections as tuples 
        (pathway_name, ligand, receptor).
    communication_clusterid : np.ndarray
        The group id of the cell-cell communications.
    D : np.ndarray
        The dissimilarity matrix for the cell-cell communications.
    

    References
    ----------
    .. [Donnat2018] Donnat, C., Zitnik, M., Hallac, D., & Leskovec, J. (2018, July). Learning structural node embeddings 
        via diffusion wavelets. In Proceedings of the 24th ACM SIGKDD International 
        Conference on Knowledge Discovery & Data Mining (pp. 1320-1329).

    """
    # Get a list of keys if only pathway_name is specified.
    if keys is None:
        keys = []
        if isinstance(pathway_name, str):
            pathways = [pathway_name]
        elif isinstance(pathway_name, list):
            pathways = pathway_name
        for pathway in pathways:
            df_ligrec = adata.uns['commot-%s-info' % pathway]['df_ligrec']
            for i in range(df_ligrec.shape[0]):
                key = (pathway, df_ligrec.iloc[i][0], df_ligrec.iloc[i][1])
                keys.append(key)
            keys.append( (pathway,'total','total') )
    nkey = len(keys)
    ncell = adata.shape[0]

    # Get a dissimilarity matrix D
    if dissimilarity_method == 'graphwave':
        if knn > 0:
            A_knn = kneighbors_graph(adata.obsm['spatial'],
                knn, mode = 'connectivity', include_self = False)
            A_spatial = A_knn + A_knn.T
            A_spatial.eliminate_zeros()
            A_spatial.data[:] = 1
        elif knn == 0:
            A_spatial = sparse.csr_matrix((ncell,ncell))
        heat_mats = []
        for key in keys:
            A_signal = adata.obsp['commot-%s-%s-%s' % (key[0], key[1], key[2])]
            A_signal_bin = binarize_sparse_matrix(A_signal, method = bin_method,
                append_zeros = bin_append_zeros, random_state = bin_random_state)
            A_signal_bin_sym = A_signal_bin + A_signal_bin.T
            A_signal_bin_sym.eliminate_zeros()
            A_signal_bin_sym.data[:] = 1
            A = A_spatial + A_signal_bin_sym
            gw = karateclub.GraphWave(**kw_graphwave)
            G = nx.from_scipy_sparse_matrix(A)
            gw.fit(G)
            R = gw.get_embedding()
            heat_mats.append(R)
        D = np.zeros([nkey,nkey],float)
        for i in range(nkey-1):
            for j in range(i+1, nkey):
                d = np.sqrt( np.linalg.norm(heat_mats[i] - heat_mats[j]) ** 2 / float(ncell) )
                D[i,j] = d; D[j,i] = d

    # Perform clustering on D
    leiden_k = min(leiden_k, len(keys)-1)
    communication_clusterid = leiden_clustering(D,
        k = leiden_k, resolution = leiden_resolution,
        random_seed = leiden_random_seed, 
        n_iterations = leiden_n_iterations)

    return keys, communication_clusterid, D


def group_communication_direction(
    adata: anndata.AnnData,
    pathway_name: Union[str, list] = None,
    keys: list = None,
    summary: str = 'sender',
    knn_smoothing: int = -1,
    normalize_vf: str = 'quantile',
    normalize_quantile: float = 0.99,
    dissimilarity_method: str = 'dot_product',
    leiden_k: int = 5,
    leiden_resolution: float = 1.0,
    leiden_random_seed: int = 1,
    leiden_n_iterations: int = -1
):
    """
    Idenfitify groups of communication directions with similar
    pattern.

    Parameters
    ----------
    adata
        The data matrix with the communication direction
        info stored in ``adata.obsm``, e.g. adata.obsm['commot_sender_vf-pathway-lig-rec'].
    pathway_name
        A string for one signaling pathway or a list of pathways.
    keys
        A list of keys for the analyzed communication connections as tuples 
        (pathway_name, ligand, receptor). 
        If given, pathway_name will be ignored.
    summary
        If 'sender', use the vector field describing to which direction the signals are sent.
        If 'receiver', use the vector field describing from which direction the signals are from.
    knn_smoothing
        The number of neighbors to smooth the communication direction. 
        If -1, no smoothing is performed.
    normalize_vf
        If 'quantile', divide all values by the length 
        given by the normalize_quantile parameter.
        If 'unit_norm', normalize each individual vector into unit norm.
        If None, original unit is used.
    normalize_quantile
        The quantile parameter to use if normalize_vf is set to 'quantile'.
    dissimilarity_method
        Currently, only dot_product is implemented.
    leiden_k
        Number of neighbors for the knn-graph to be fed to leiden clustering algorithm.
    leiden_resolution
        The resolution parameter for the leiden clustering algorithm.
    leiden_random_seed
        The random seed for the leiden clustering algorithm.
    leiden_n_iterations
        The maximum number of iterations for the leiden algorithm.
        The algorithm will run until convergence if set to -1.

    Returns
    -------
    keys : list
        The list of keys for the analyzed communication connections as tuples 
        (pathway_name, ligand, receptor).
    direction_clusterid : np.ndarray
        The group id of the communication directions.
    D : np.ndarray
        The dissimilarity matrix for the communication directions.
    
    """
    # Get a list of keys if only pathway is given.
    if keys is None:
        keys = []
        if isinstance(pathway_name, str):
            pathways = [pathway_name]
        elif isinstance(pathway_name, list):
            pathways = pathway_name
        for pathway in pathways:
            df_ligrec = adata.uns['commot-%s-info' % pathway]['df_ligrec']
            for i in range(df_ligrec.shape[0]):
                key = (pathway, df_ligrec.iloc[i][0], df_ligrec.iloc[i][1])
                keys.append(key)
            keys.append( (pathway,'total','total') )
    # Process the vector fields
    V_list = []
    for key in keys:
        V = adata.obsm['commot_%s_vf-%s-%s-%s' % (summary, key[0], key[1], key[2])]
        V_processed = preprocess_vector_field(adata.obsm['spatial'],
            V, knn_smoothing = knn_smoothing, normalize_vf = normalize_vf,
            quantile = normalize_quantile)
        V_list.append(V_processed)
    # Get a distance matrix between the vector fields
    nV = len(keys)
    D = np.zeros([nV,nV], float)
    for i in range(nV-1):
        Vi = V_list[i]
        for j in range(i+1,nV):
            Vj = V_list[j]
            if dissimilarity_method == 'dot_product':
                d = np.exp( - ( ( Vi * Vj ).sum(axis=1) ).mean() )
            D[i,j] = d; D[j,i] = d
    # Cluster the vector fields with the D matrix
    # Perform clustering
    leiden_k = min(leiden_k, len(keys)-1)
    direction_clusterid = leiden_clustering(D,
        k = leiden_k, resolution = leiden_resolution,
        random_seed = leiden_random_seed, 
        n_iterations = leiden_n_iterations)
    
    return keys, direction_clusterid, D


def communication_spatial_autocorrelation(
    adata: anndata.AnnData,
    pathway_name: Union[str, list] = None,
    keys: list = None,
    method: str = 'Moran',
    normalize_vf: bool = False,
    summary: str = 'sender',
    weight_bandwidth: float = None,
    weight_k: int = 10,
    weight_function: str = 'triangular',
    weight_row_standardize: bool = False,
    n_permutations: int = 999
):
    """
    Spatial autocorrelation of communication directions.

    Parameters
    ----------
    adata
        The data matrix with the communication vector fields
        info stored in ``adata.ubsm``.
    pathway_name
        A string for one signaling pathway or a list of pathways.
    keys
        A list of keys for the analyzed communication connections as tuples 
        (pathway_name, ligand, receptor). 
        If given, pathway_name will be ignored.
    method
        The method to use. Currently, only Moran's I [Liu2015]_ for vectors is implemented.
    normalize_vf
        Whether to normalize the vector field so that the autocorrelation only reflects
        directions.
    summary
        If 'sender', use the vector field describing to which direction the signals are sent.
        If 'receiver', use the vector field describing from which direction the signals are from.
    weight_bandwidth
        The bandwidth for the kernel to assign knn graph weights.
        If given, weight_k is ignored.
    weight_k
        The number of nearest neighbors for the knn graph.
    weight_function
        Kernel functions for assigning weight. 
        Choices: 'triangular','uniform','quadratic','quartic','gaussian'.
        See libpysal.weights.Kernel of the ``libpysal`` package for details.
    weight_row_standardize
        Whether to standardize the weights so that the heterogeneity in local cell/position
        density does not affect the results.
    n_permutations
        Number of permutations for computing p-values.

    Returns
    -------
    keys : list
        The list of keys for the analyzed communication connections as tuples 
        (pathway_name, ligand, receptor).
    moranI : np.ndarray
        A vector of moran's I statistics for corresponding to each key in keys.
    p_value : np.ndarray
        The p-values.

    References
    ----------
    .. [Liu2015] Liu, Y., Tong, D., & Liu, X. (2015). Measuring spatial 
        autocorrelation of vectors. Geographical Analysis, 47(3), 300-319.

    """
    if keys is None:
        keys = []
        if isinstance(pathway_name, str):
            pathways = [pathway_name]
        elif isinstance(pathway_name, list):
            pathways = pathway_name
        for pathway in pathways:
            df_ligrec = adata.uns['commot-%s-info' % pathway]['df_ligrec']
            for i in range(df_ligrec.shape[0]):
                key = (pathway, df_ligrec.iloc[i][0], df_ligrec.iloc[i][1])
                keys.append(key)
            keys.append( (pathway,'total','total') )
    
    moranI = []
    p_value = []
    X = adata.obsm['spatial']
    for key in keys:
        V = adata.obsm['commot_%s_vf-%s-%s-%s' % (summary, key[0], key[1], key[2])]
        if normalize_vf:
            V = normalize(V)
        I,p = moranI_vector_global(X, V,
            weight_bandwidth = weight_bandwidth,
            weight_k = weight_k,
            weight_function = weight_function,
            weight_row_standardize = weight_row_standardize,
            n_permutations = n_permutations)
        moranI.append(I)
        p_value.append(p)
    moranI = np.array(moranI, float)
    p_value = np.array(p_value, float)

    return keys, moranI, p_value


