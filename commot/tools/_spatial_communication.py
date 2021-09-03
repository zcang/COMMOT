from typing import Optional
import gc
import ot
import sys
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import cm
import matplotlib.pyplot as plt
import plotly
from scipy import sparse
from scipy.spatial import distance_matrix
from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import normalize

from .._optimal_transport import cot_sparse
from .._optimal_transport import cot_combine_sparse
from .._optimal_transport import uot
from .._optimal_transport import usot

def kernel_function(x, eta, nu, kernel, normalization=None):
    if kernel == "exp":
        phi = np.exp(-np.power(x/eta, nu))
    elif kernel == "lorentz":
        phi = 1 / ( 1 + np.power(x/eta, nu) )
    if normalization == "unit_row_sum":
        phi = (phi.T / np.sum(phi.T, axis=0)).T
    elif normalization == "unit_col_sum":
        phi = phi / np.sum(phi, axis=0)
    return phi

def coo_from_dense_submat(row, col, x, shape):
    col_ind, row_ind = np.meshgrid(col, row)
    sparse_mat = sparse.coo_matrix((x.reshape(-1), (row_ind.reshape(-1), col_ind.reshape(-1))), shape=shape)
    sparse_mat.eliminate_zeros()
    return(sparse_mat)

class CellCommunication(object):

    def __init__(self,
        adata,
        df_ligrec,
        dmat,
        dis_thr,
        cost_scale,
        cost_type
    ):
        data_genes = set(adata.var_names)
        self.ligs = list(set(df_ligrec[0]).intersection(data_genes))
        self.recs = list(set(df_ligrec[1]).intersection(data_genes))
        A = np.inf * np.ones([len(self.ligs), len(self.recs)], float)
        for i in range(len(df_ligrec)):
            tmp_lig, tmp_rec = df_ligrec.iloc[i]
            if tmp_lig in self.ligs and tmp_rec in self.recs:
                if cost_scale is None:
                    A[self.ligs.index(tmp_lig), self.recs.index(tmp_rec)] = 1.0
                else:
                    A[self.ligs.index(tmp_lig), self.recs.index(tmp_rec)] = cost_scale[(tmp_lig, tmp_rec)]
        self.A = A
        self.S = adata[:,self.ligs].X.toarray()
        self.D = adata[:,self.recs].X.toarray()
        if cost_type == 'euc':
            self.M = dmat
        elif cost_type == 'euc_square':
            self.M = dmat ** 2
        if np.isscalar(dis_thr):
            if cost_type == 'euc_square':
                dis_thr = dis_thr ** 2
            self.cutoff = float(dis_thr) * np.ones_like(A)
        elif type(dis_thr) is dict:
            self.cutoff = np.zeros_like(A)
            for i in range(A.shape[0]):
                for j in range(A.shape[1]):
                    if A[i,j] > 0:
                        if cost_type == 'euc_square':
                            self.cutoff[i,j] = dis_thr[(self.ligs[i], self.recs[j])] ** 2
                        else:
                            self.cutoff[i,j] = dis_thr[(self.ligs[i], self.recs[j])]
        self.nlig = self.S.shape[1]; self.nrec = self.D.shape[1]
        self.npts = adata.shape[0]

    def run_cot_signaling(self,
        cot_eps_p=1e-1, 
        cot_eps_mu=None, 
        cot_eps_nu=None, 
        cot_rho=1e1, 
        cot_nitermax=1e4, 
        cot_weights=(0.25,0.25,0.25,0.25), 
        smooth=False, 
        smth_eta=None, 
        smth_nu=None, 
        smth_kernel=None
    ):
        if not smooth:
            self.comm_network = cot_combine_sparse(self.S, self.D, self.A, self.M, self.cutoff, \
                eps_p=cot_eps_p, eps_mu=cot_eps_mu, eps_nu=cot_eps_nu, rho=cot_rho, weights=cot_weights, nitermax=cot_nitermax)
        if smooth:
            # Smooth the expression
            S_smth = np.zeros_like(self.S)
            D_smth = np.zeros_like(self.D)
            for i in range(self.nlig):
                nzind = np.where(self.S[:,i] > 0)[0]
                phi = kernel_function(self.M[nzind,:][:,nzind], smth_eta, smth_nu, smth_kernel, normalization='unit_col_sum')
                S_smth[nzind,i] = np.matmul( phi, self.S[nzind,i].reshape(-1,1) )[:,0]
            for i in range(self.nrec):
                nzind = np.where(self.D[:,i] > 0)[0]
                phi = kernel_function(self.M[nzind,:][:,nzind], smth_eta, smth_nu, smth_kernel, normalization='unit_col_sum')
                D_smth[nzind,i] = np.matmul( phi, self.D[nzind,i].reshape(-1,1) )[:,0]
            # Get the transport plans for the smoothed distributions
            P_smth = cot_combine_sparse(S_smth, D_smth, self.A, self.M, self.cutoff, \
                eps_p=cot_eps_p, eps_mu=cot_eps_mu, eps_nu=cot_eps_nu, rho=cot_rho, weights=cot_weights, nitermax=cot_nitermax)
            # Deconvolute back to original cells
            self.comm_network = {}
            for i in range(self.nlig):
                S = self.S[:,i]
                nzind_S = np.where(S > 0)[0]
                phi_S = kernel_function(self.M[nzind_S,:][:,nzind_S], smth_eta, smth_nu, smth_kernel, normalization='unit_col_sum')
                S_contrib = phi_S * S[nzind_S]; S_contrib = S_contrib / np.sum(S_contrib, axis=1, keepdims=True)
                for j in range(self.nrec):
                    D = self.D[:,j]
                    nzind_D = np.where(D > 0)[0]
                    if np.isinf(self.A[i,j]): continue
                    P_dense = P_smth[(i,j)].toarray()
                    P_sub = P_dense[nzind_S,:][:,nzind_D]
                    phi_D = kernel_function(self.M[nzind_D,:][:,nzind_D], smth_eta, smth_nu, smth_kernel, normalization='unit_col_sum')
                    D_contrib = phi_D * D[nzind_D]; D_contrib = D_contrib / np.sum(D_contrib, axis=1, keepdims=True)
                    P_deconv = np.matmul(S_contrib.T, np.matmul(P_sub, D_contrib))
                    for k in range(len(nzind_D)):
                        P_dense[nzind_S, nzind_D[k]] = P_deconv[:,k]
                    P_deconv_sparse = coo_from_dense_submat(nzind_S, nzind_D, P_deconv, shape=(self.npts, self.npts))
                    self.comm_network[(i,j)] = P_deconv_sparse

    def smooth(self, eta, nu, kernel):
        S_smth = np.zeros_like(self.S)
        D_smth = np.zeros_like(self.D)
        for i in range(self.nlig):
            nzind = np.where(self.S[:,i] > 0)[0]
            phi = kernel_function(self.M[nzind,:][:,nzind], eta, nu, kernel, normalization='unit_col_sum')
            S_smth[nzind,i] = np.matmul( phi, self.S[nzind,i].reshape(-1,1) )[:,0]
        for i in range(self.nrec):
            nzind = np.where(self.D[:,i] > 0)[0]
            phi = kernel_function(self.M[nzind,:][:,nzind], eta, nu, kernel, normalization='unit_col_sum')
            D_smth[nzind,i] = np.matmul( phi, self.D[nzind,i].reshape(-1,1) )[:,0]
        self.S_smth = S_smth
        self.D_smth = D_smth
        self.kernel_eta = eta
        self.kernel_nu = nu
        self.kernel = kernel

def assign_distance(adata, dmat=None):
    if dmat is None:
        adata.obsp["spatial_distance"] = distance_matrix(adata.obsm["spatial"], adata.obsm["spatial"])
    else:
        adata.obsp["spatial_distance"] = dmat

def summarize_cluster(X, clusterid, clusternames, n_permutations=500):
    # Input a sparse matrix of cell signaling and output a pandas dataframe
    # for cluster-cluster signaling
    n = len(clusternames)
    X_cluster = np.empty([n,n], float)
    p_cluster = np.zeros([n,n], float)
    for i in range(n):
        tmp_idx_i = np.where(clusterid==clusternames[i])[0]
        for j in range(n):
            tmp_idx_j = np.where(clusterid==clusternames[j])[0]
            X_cluster[i,j] = X[tmp_idx_i,:][:,tmp_idx_j].mean()
    for i in range(n_permutations):
        clusterid_perm = np.random.permutation(clusterid)
        X_cluster_perm = np.empty([n,n], float)
        for j in range(n):
            tmp_idx_j = np.where(clusterid_perm==clusternames[j])[0]
            for k in range(n):
                tmp_idx_k = np.where(clusterid_perm==clusternames[k])[0]
                X_cluster_perm[j,k] = X[tmp_idx_j,:][:,tmp_idx_k].mean()
        p_cluster[X_cluster_perm >= X_cluster] += 1.0
    p_cluster = p_cluster / n_permutations
    df_cluster = pd.DataFrame(data=X_cluster, index=clusternames, columns=clusternames)
    df_p_value = pd.DataFrame(data=p_cluster, index=clusternames, columns=clusternames)
    return df_cluster, df_p_value

def cluster_center(adata, clustering, method="geometric_mean"):
    X = adata.obsm['spatial']
    cluster_pos = {}
    clusternames = list( adata.obs[clustering].unique() )
    clusternames.sort()
    clusterid = np.array( adata.obs[clustering], str )
    for name in clusternames:
        tmp_idx = np.where(clusterid==name)[0]
        tmp_X = X[tmp_idx]
        if method == "geometric_mean":
            X_mean = np.mean(tmp_X, axis=0)
        elif method == "representative_point":
            tmp_D = distance_matrix(tmp_X, tmp_X)
            tmp_D = tmp_D ** 2
            X_mean = tmp_X[np.argmin(tmp_D.sum(axis=1)),:]
        cluster_pos[name] = X_mean
    return cluster_pos

def spatial_communication(
    adata: anndata.AnnData, 
    pathway_name: str = None, 
    df_ligrec: pd.DataFrame = None, 
    dis_thr: Optional[float] = None, 
    cost_scale: Optional[dict] = None, 
    cost_type: str = 'euc',
    cot_eps_p: float = 1e-1, 
    cot_eps_mu: Optional[float] = None, 
    cot_eps_nu: Optional[float] = None, 
    cot_rho: float =1e1, 
    cot_nitermax: int = 10000, 
    cot_weights: tuple = (0.25,0.25,0.25,0.25), 
    smooth: bool = False, 
    smth_eta: float = None, 
    smth_nu: float = None, 
    smth_kernel: str = 'exp',
    copy: bool = False
):
    """
    Infer spatial communication.

    Parameters
    ----------
    adata
        The data matrix of shape ``n_obs`` × ``n_var``.
        Rows correspond to cells or positions and columns to genes.
        If the spatial distance is absent in ``.obsp['spatial_distance']``, Euclidean distance determined from ``.obsm['spatial']`` will be used.
    pathway_name
        Name of the signaling pathway.
    df_ligrec
        A data frame where each row corresponds to a ligand-receptor pair with ligands and receptors in the first and second columns, respectively.
    dis_thr
        The threshold of spatial distance of signaling.
    cost_scale
        Weight coefficients of the cost matrix for each ligand-receptor pair, e.g. cost_scale[('ligA','recA')] specifies weight for the pair ligA and recA.
        If None, all pairs have the same weight. 
    cost_type
        If 'euc', the original Euclidean distance will be used as cost matrix. If 'euc_square', the square of the Euclidean distance will be used.
    cot_eps_p
        The coefficient of entropy regularization for transport plan.
    cot_eps_mu
        The coefficient of entropy regularization for untransported source (ligand). Set to equal to cot_eps_p for fast algorithm.
    cot_eps_nu
        The coefficient of entropy regularization for unfulfilled target (receptor). Set to equal to cot_eps_p for fast algorithm.
    cot_rho
        The coefficient of penalty for unmatched mass.
    cot_nitermax
        Maximum iteration for collective optimal transport algorithm.
    cot_weights
        A tuple of four weights that add up to one. The weights corresponds to four setups of collective optimal transport: 
        1) all ligands-all receptors, 2) each ligand-all receptors, 3) all ligands-each receptor, 4) each ligand-each receptor.
    smooth
        Whether to (spatially) smooth the gene expression for identifying more global signaling trend.
    smth_eta
        Kernel bandwidth for smoothing
    smth_nu
        Kernel sharpness for smoothing
    smth_kernel
        'exp' exponential kernel. 'lorentz' Lorentz kernel.
    copy
        Whether to return a copy of the :class:`anndata.AnnData`.

    Returns
    -------
    adata : anndata.AnnData
        Signaling matrices are added to ``.obsp``, e.g., for a signaling pathway named "pathway", 
        ``.obsp['commot-pathway-ligA-recA']``
        is a ``n_obs`` × ``n_obs`` matrix with the *ij* th entry being the "score" of 
        cell *i* sending signal to cell *j* through ligA and recA.
        The marginal sums of the signaling matrices are stored in ``.obsm['commot-pathway-sum']``.
        Metadata of the analysis is added to ``.uns['commot-pathway-info']``.
        If copy=True, return the AnnData object and return None otherwise.

    """

    if not 'spatial_distance' in adata.obsp.keys():
        dis_mat = distance_matrix(adata.obsm["spatial"], adata.obsm["spatial"])
        # assign_distance(adata, dmat=None)
    else:
        dis_mat = adata.obsp['spatial_distance']

    # remove unavailable genes from df_ligrec
    data_genes = list(adata.var_names)
    tmp_ligrec = []
    for i in range(df_ligrec.shape[0]):
        if df_ligrec.iloc[i][0] in data_genes and df_ligrec.iloc[i][1] in data_genes:
            tmp_ligrec.append([df_ligrec.iloc[i][0], df_ligrec.iloc[i][1]])
    tmp_ligrec = np.array(tmp_ligrec, str)
    df_ligrec = pd.DataFrame(data=tmp_ligrec)
    # Drop duplicate pairs
    df_ligrec = df_ligrec.drop_duplicates()


    model = CellCommunication(adata, 
        df_ligrec, 
        dis_mat, 
        dis_thr, 
        cost_scale, 
        cost_type
    )
    model.run_cot_signaling(cot_eps_p=1e-1, 
        cot_eps_mu = cot_eps_mu, 
        cot_eps_nu = cot_eps_nu, 
        cot_rho = cot_rho, 
        cot_nitermax = cot_nitermax, 
        cot_weights = cot_weights, 
        smooth = smooth, 
        smth_eta = smth_eta, 
        smth_nu = smth_nu, 
        smth_kernel = smth_kernel
    )

    
    adata.uns['commot-'+pathway_name+'-info'] = {}
    df_ligrec_write = pd.DataFrame(data=df_ligrec.values, columns=['ligand','receptor'])
    adata.uns['commot-'+pathway_name+'-info']['df_ligrec'] = df_ligrec_write
    adata.uns['commot-'+pathway_name+'-info']['distance_threshold'] = dis_thr

    ncell = adata.shape[0]
    X_sender = np.empty([ncell,0], float)
    X_receiver = np.empty([ncell,0], float)
    col_names_sender = []
    col_names_receiver = []
    tmp_ligs = model.ligs
    tmp_recs = model.recs
    S_total = sparse.csr_matrix((ncell, ncell), dtype=float)
    for (i,j) in model.comm_network.keys():
        S = model.comm_network[(i,j)]
        adata.obsp['commot-'+pathway_name+'-'+tmp_ligs[i]+'-'+tmp_recs[j]] = S
        S_total = S_total + S
        lig_sum = np.array(S.sum(axis=1))
        rec_sum = np.array(S.sum(axis=0).T)
        X_sender = np.concatenate((X_sender, lig_sum), axis=1)
        X_receiver = np.concatenate((X_receiver, rec_sum), axis=1)
        col_names_sender.append("sender-%s-%s" % (tmp_ligs[i], tmp_recs[j]))
        col_names_receiver.append("receiver-%s-%s" % (tmp_ligs[i], tmp_recs[j]))
    X_sender = np.concatenate((X_sender, X_sender.sum(axis=1).reshape(-1,1)), axis=1)
    X_receiver = np.concatenate((X_receiver, X_receiver.sum(axis=1).reshape(-1,1)), axis=1)
    col_names_sender.append("sender-total-total")
    col_names_receiver.append("receiver-total-total")
    X = np.concatenate((X_sender, X_receiver), axis=1)
    col_names = col_names_sender + col_names_receiver    

    adata.obsp['commot-'+pathway_name+'-total-total'] = S_total

    df = pd.DataFrame(data=X, columns=col_names, index=adata.obs_names)
    adata.obsm['commot-'+pathway_name+"-sum"] = df

    del model
    gc.collect()

    return adata if copy else None


def communication_direction(
    adata: anndata.AnnData,
    pathway_name: str = None,
    lr_pair = None,
    k: int = 5,
    pos_idx: Optional[np.ndarray] = None,
    copy: bool = False
):
    """
    Construct spatial vector fields for inferred communication.

    Parameters
    ----------
    adata
        The data matrix of shape ``n_obs`` × ``n_var``.
        Rows correspond to cells or positions and columns to genes.
    pathway_name
        Name of the signaling pathway.
    lr_pair
        A tuple of ligand-receptor pair. 
        If given, only the signaling direction of this pair is computed.
    k
        Top k senders or receivers to consider when determining the direction.
    pos_idx
        The columns in ``.obsm['spatial']`` to use. If None, all columns are used.
    copy
        Whether to return a copy of the :class:`anndata.AnnData`.
    
    Returns
    -------
    adata : anndata.AnnData
        Vector fields describing signaling directions are added to ``.obsm``, 
        e.g., for a signaling pathway named "pathway", 
        ``.obsm['commot_sender_vf-pathway-ligA-recA']`` and ``.obsm['commot_receiver_vf-pathway-ligA-recA']``
        describe the signaling directions of the cells as, respectively, senders and receivers through the 
        ligand-receptor pair ligA and recA.
        If copy=True, return the AnnData object and return None otherwise.

    """

    name_mat = adata.uns['commot-'+pathway_name+'-info']['df_ligrec'].values
    name_mat = np.concatenate((name_mat, np.array([['total','total']],str)), axis=0)
    if not lr_pair is None:
        name_mat = np.array([[lr_pair[0], lr_pair[1]]], dtype=str)
    ncell = adata.shape[0]
    pts = np.array( adata.obsm['spatial'], float )
    if not pos_idx is None:
        pts = pts[:,pos_idx]
    storage = 'sparse'

    if storage == 'dense':
        for i in range(name_mat.shape[0]):
            lig = name_mat[i,0]
            rec = name_mat[i,1]
            S_np = adata.obsp['commot-'+pathway_name+'-'+lig+'-'+rec].toarray()
            sender_vf = np.zeros_like(pts)
            receiver_vf = np.zeros_like(pts)

            tmp_idx = np.argsort(-S_np,axis=1)[:,:k]
            avg_v = np.zeros_like(pts)
            for ik in range(k):
                tmp_v = pts[tmp_idx[:,ik]] - pts[np.arange(ncell,dtype=int)]
                tmp_v = normalize(tmp_v, norm='l2')
                avg_v = avg_v + tmp_v * S_np[np.arange(ncell,dtype=int),tmp_idx[:,ik]].reshape(-1,1)
            avg_v = normalize(avg_v)
            sender_vf = avg_v * np.sum(S_np,axis=1).reshape(-1,1)

            S_np = S_np.T
            tmp_idx = np.argsort(-S_np,axis=1)[:,:k]
            avg_v = np.zeros_like(pts)
            for ik in range(k):
                tmp_v = -pts[tmp_idx[:,ik]] + pts[np.arange(ncell,dtype=int)]
                tmp_v = normalize(tmp_v, norm='l2')
                avg_v = avg_v + tmp_v * S_np[np.arange(ncell,dtype=int),tmp_idx[:,ik]].reshape(-1,1)
            avg_v = normalize(avg_v)
            receiver_vf = avg_v * np.sum(S_np,axis=1).reshape(-1,1)

            del S_np

    elif storage == 'sparse':
        for i in range(name_mat.shape[0]):
            lig = name_mat[i,0]
            rec = name_mat[i,1]
            S = adata.obsp['commot-'+pathway_name+'-'+lig+'-'+rec]
            S_sum_sender = np.array( S.sum(axis=1) ).reshape(-1)
            S_sum_receiver = np.array( S.sum(axis=0) ).reshape(-1)
            sender_vf = np.zeros_like(pts)
            receiver_vf = np.zeros_like(pts)

            S_lil = S.tolil()
            for j in range(S.shape[0]):
                if len(S_lil.rows[j]) <= k:
                    tmp_idx = np.array( S_lil.rows[j], int )
                    tmp_data = np.array( S_lil.data[j], float )
                else:
                    row_np = np.array( S_lil.rows[j], int )
                    data_np = np.array( S_lil.data[j], float )
                    sorted_idx = np.argsort( -data_np )[:k]
                    tmp_idx = row_np[ sorted_idx ]
                    tmp_data = data_np[ sorted_idx ]
                if len(tmp_idx) == 0:
                    continue
                elif len(tmp_idx) == 1:
                    avg_v = pts[tmp_idx[0],:] - pts[j,:]
                else:
                    tmp_v = pts[tmp_idx,:] - pts[j,:]
                    tmp_v = normalize(tmp_v, norm='l2')
                    avg_v = tmp_v * tmp_data.reshape(-1,1)
                    avg_v = np.sum( avg_v, axis=0 )
                avg_v = normalize( avg_v.reshape(1,-1) )
                sender_vf[j,:] = avg_v[0,:] * S_sum_sender[j]
            
            S_lil = S.T.tolil()
            for j in range(S.shape[0]):
                if len(S_lil.rows[j]) <= k:
                    tmp_idx = np.array( S_lil.rows[j], int )
                    tmp_data = np.array( S_lil.data[j], float )
                else:
                    row_np = np.array( S_lil.rows[j], int )
                    data_np = np.array( S_lil.data[j], float )
                    sorted_idx = np.argsort( -data_np )[:k]
                    tmp_idx = row_np[ sorted_idx ]
                    tmp_data = data_np[ sorted_idx ]
                if len(tmp_idx) == 0:
                    continue
                elif len(tmp_idx) == 1:
                    avg_v = -pts[tmp_idx,:] + pts[j,:]
                else:
                    tmp_v = -pts[tmp_idx,:] + pts[j,:]
                    tmp_v = normalize(tmp_v, norm='l2')
                    avg_v = tmp_v * tmp_data.reshape(-1,1)
                    avg_v = np.sum( avg_v, axis=0 )
                avg_v = normalize( avg_v.reshape(1,-1) )
                receiver_vf[j,:] = avg_v[0,:] * S_sum_receiver[j]



            adata.obsm["commot_sender_vf-"+pathway_name+'-'+lig+'-'+rec] = sender_vf
            adata.obsm["commot_receiver_vf-"+pathway_name+'-'+lig+'-'+rec] = receiver_vf

    return adata if copy else None

def cluster_communication(
    adata: anndata.AnnData,
    pathway_name: str = None,
    clustering: str = None,
    n_permutations: int = 100,
    copy: bool = False
):
    """
    Summarize cell-cell communication to cluster-cluster communication.

    Parameters
    ----------
    adata
        The data matrix of shape ``n_obs`` × ``n_var``.
        Rows correspond to cells or positions and columns to genes.
    pathway_name
        Name of the signaling pathway.
    clustering
        Name of clustering with the labels stored in ``.obs[clustering]``.
    n_permutations
        Number of label permutations for computing the p-value.
    copy
        Whether to return a copy of the :class:`anndata.AnnData`.
    
    Returns
    -------
    adata : anndata.AnnData
        Add cluster-cluster communication matrix to 
        ``.uns['commot_cluster-pathway_name-clustering-ligA-recA']``
        for the signaling pathway named 'pathway_name' and the cell clustering 
        named 'clustering' through the ligand-receptor pair 'ligA' and 'recA'.
        The first object is the communication matrix and the second object contains
        the corresponding p-values.
        If copy=True, return the AnnData object and return None otherwise.

    """
    celltypes = list( adata.obs[clustering].unique() )
    celltypes.sort()
    clusterid = np.array(adata.obs[clustering], str)
    name_mat = adata.uns['commot-'+pathway_name+'-info']['df_ligrec'].values
    name_mat = np.concatenate((name_mat, np.array([['total','total']],str)), axis=0)
    for i in range(name_mat.shape[0]):
        lig = name_mat[i,0]
        rec = name_mat[i,1]
        S = adata.obsp['commot-'+pathway_name+'-'+lig+'-'+rec]
        tmp_df, tmp_p_value = summarize_cluster(S,
            clusterid, celltypes, n_permutations=n_permutations)
        adata.uns['commot_cluster-'+pathway_name+'-'+clustering+'-'+lig+'-'+rec] = {'communication_matrix': tmp_df, 'communication_pvalue': tmp_p_value}
    
    return adata if copy else None

def cluster_position(
    adata: anndata.AnnData,
    clustering: str = None,
    method: str = 'geometric_mean',
    copy: bool = False
):
    """
    Assign spatial positions to clusters.

    Parameters
    ----------
    adata
        The data matrix of shape ``n_obs`` × ``n_var``.
        Rows correspond to cells or positions and columns to genes.
    clustering
        Name of clustering with the labels stored in ``.obs[clustering]``.
    method
        'geometric_mean' geometric mean of the positions. 
        'representative_point' the position in the cluster with 
        minimum total distance to other points.
    copy
        Whether to return a copy of the :class:`anndata.AnnData`.

    Returns
    -------
    adata : anndata.AnnData
        Add cluster positions to ``.uns['cluster_pos-clustering_name']`` for the clustering named
        'clustering_name'.
        If copy=True, return the AnnData object and return None otherwise.

    """
    cluster_pos = cluster_center(adata, clustering, method=method)
    adata.uns['cluster_pos-'+clustering] = cluster_pos

    return adata if copy else None