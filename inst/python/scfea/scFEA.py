# -*- coding: utf-8 -*-
"""
scFEA: single-cell Flux Estimation Analysis.

Vendored into scop from https://github.com/changwn/scFEA
(commit 4c1fb76d52f07bafad84ce7686ad7c3acfcf0126).

Original author: wnchang@iu.edu
Reference: Alghamdi et al. (2021), Genome Research.

This module keeps the original GNN architecture, loss and training loop
bit-faithful to the upstream ``main()``. Only the I/O has been refactored:
the expression matrix is accepted as an in-memory pandas DataFrame and the
results are returned as DataFrames instead of being read from / written to
CSV files. The model files are loaded from a ``data_dir`` prepared by the R
wrapper, keeping data resources outside the R package tarball.
"""

# system lib
import os
import time
import warnings

# tools
import torch
from torch.autograd import Variable
import numpy as np
import pandas as pd
from tqdm import tqdm

# scFEA lib
from .ClassFlux import FLUX  # Flux class network
from .util import pearsonr
from .DatasetFlux import MyDataset


# hyper parameters
LEARN_RATE = 0.008
LAMB_BA = 1
LAMB_NG = 1
LAMB_CELL = 1
LAMB_MOD = 1e-2


# bundled M168 asset sets, keyed by species
_SPECIES_FILES = {
    'human': {
        'moduleGene_file': 'module_gene_m168.csv',
        'stoichiometry_matrix': 'cmMat_c70_m168.csv',
        'cName_file': 'cName_c70_m168.csv',
    },
    'mouse': {
        'moduleGene_file': 'module_gene_complete_mouse_m168.csv',
        'stoichiometry_matrix': 'cmMat_complete_mouse_c70_m168.csv',
        'cName_file': 'cName_complete_mouse_c70_m168.csv',
    },
}


def myLoss(m, c, lamb1=0.2, lamb2=0.2, lamb3=0.2, lamb4=0.2,
           geneScale=None, moduleScale=None):

    # balance constrain
    total1 = torch.pow(c, 2)
    total1 = torch.sum(total1, dim=1)

    # non-negative constrain
    error = torch.abs(m) - m
    total2 = torch.sum(error, dim=1)

    # sample-wise variation constrain
    diff = torch.pow(torch.sum(m, dim=1) - geneScale, 2)
    if sum(diff > 0) == m.shape[0]:  # solve Nan after several iteraions
        total3 = torch.pow(diff, 0.5)
    else:
        total3 = diff

    # module-wise variation constrain
    if lamb4 > 0:
        corr = torch.ones(m.shape[0], dtype=m.dtype, device=m.device)
        for i in range(m.shape[0]):
            corr[i] = pearsonr(m[i, :], moduleScale[i, :])
        corr = torch.abs(corr)
        penal_m_var = torch.ones(m.shape[0], dtype=m.dtype, device=m.device) - corr
        total4 = penal_m_var
    else:
        total4 = torch.zeros(m.shape[0], dtype=m.dtype, device=m.device)

    # loss
    loss1 = torch.sum(lamb1 * total1)
    loss2 = torch.sum(lamb2 * total2)
    loss3 = torch.sum(lamb3 * total3)
    loss4 = torch.sum(lamb4 * total4)
    loss = loss1 + loss2 + loss3 + loss4
    return loss, loss1, loss2, loss3, loss4


def _build_X_batch(geneExpr_np, module_gene_indices, n_genes, n_modules,
                   cell_indices, device):
    """Build the expanded X tensor for a subset of cells, directly as a torch
    tensor without intermediate pandas DataFrames.

    Parameters
    ----------
    geneExpr_np : np.ndarray (n_cells, n_genes), float32
    module_gene_indices : list of list of int
        module_gene_indices[i] = [global_gene_idx, ...] for module i.
    n_genes : int
    n_modules : int
    cell_indices : np.ndarray or slice
        Which cells (rows) to include.
    device : torch.device

    Returns
    -------
    X : torch.FloatTensor (len(cell_indices), n_modules * n_genes)
    """
    n_cells_batch = len(cell_indices)
    X = torch.zeros(n_cells_batch, n_modules * n_genes,
                    dtype=torch.float32, device=device)
    expr_slice = geneExpr_np[cell_indices, :]  # (batch, n_genes)
    for i in range(n_modules):
        indices = module_gene_indices[i]
        if not indices:
            continue
        for local_j, global_j in enumerate(indices):
            col = i * n_genes + global_j
            X[:, col] = torch.from_numpy(expr_slice[:, global_j]).float().to(device)
    return X


def _build_module_scale(geneExpr_np, module_gene_indices, moduleLen,
                        cell_indices, device):
    """Compute per-cell mean module expression for a subset of cells.

    Parameters
    ----------
    geneExpr_np : np.ndarray (n_cells, n_genes), float32
    module_gene_indices : list of list of int
    moduleLen : np.ndarray (n_modules,)  number of genes per module
    cell_indices : np.ndarray or slice
    device : torch.device

    Returns
    -------
    module_scale : torch.FloatTensor (len(cell_indices), n_modules)
    """
    n_cells_batch = len(cell_indices)
    n_modules = len(module_gene_indices)
    module_scale = torch.zeros(n_cells_batch, n_modules,
                               dtype=torch.float32, device=device)
    expr_slice = geneExpr_np[cell_indices, :]
    for i in range(n_modules):
        indices = module_gene_indices[i]
        if not indices or moduleLen[i] == 0:
            continue
        module_scale[:, i] = torch.from_numpy(
            expr_slice[:, indices].sum(axis=1) / moduleLen[i]
        ).float().to(device)
    return module_scale


def run_scfea(expr, *, species='human', sc_imputation=False, n_epoch=100,
              seed=16, verbose=True, device=None, data_dir=None,
              max_cells=None, predict_batch_size=10000):
    """Estimate single-cell metabolic flux with scFEA.

    Parameters
    ----------
    expr : pandas.DataFrame
        Gene-expression matrix, cells x genes (rows = cells, columns = gene
        symbols). May be raw or normalised counts; a log2 transform is applied
        automatically when the maximum value exceeds 50, mirroring upstream.
    species : {'human', 'mouse'}, default 'human'
        Selects the bundled M168 module / stoichiometry / compound-name files.
    sc_imputation : bool, default False
        Whether to MAGIC-impute the expression matrix before flux estimation
        (recommended for sparse 10x data). Requires the ``magic-impute``
        package; it is imported only when this is True.
    n_epoch : int, default 100
        Number of training epochs for the flux GNN. Must be > 0.
    seed : int, default 16
        Torch random seed. The default matches the upstream scFEA script.
    verbose : bool, default True
        Print progress messages and show the training progress bar.
    device : str or torch.device, optional
        Compute device. Defaults to ``cuda:0`` if available, else ``cpu``.
    data_dir : str
        Directory containing scFEA M168 resource CSV files.
    max_cells : int, optional
        Maximum number of cells used for GNN training. When the input has more
        cells, a random subset is sampled for training and the trained model
        is used to predict fluxes for all cells in batches. This dramatically
        reduces peak memory for large datasets. Default (None) trains on all
        cells, matching the original upstream behaviour.
    predict_batch_size : int, default 10000
        Number of cells processed per forward pass during the prediction
        (inference) phase. Lower this if you encounter out-of-memory errors
        during prediction.

    Returns
    -------
    dict
        ``{'flux': DataFrame (cells x modules),
        'balance': DataFrame (cells x metabolites),
        'predictions': dict}`` where ``flux`` corresponds to scFEA's ``setF``
        output, ``balance`` to its ``setB`` output, and ``predictions`` bundles
        the raw numpy arrays plus cell / module / compound name indices.
    """
    if n_epoch <= 0:
        raise NameError('n_epoch must greater than 1!')

    if species not in _SPECIES_FILES:
        raise ValueError(
            "species must be one of %s" % sorted(_SPECIES_FILES))

    if not isinstance(expr, pd.DataFrame):
        raise TypeError('expr must be a pandas DataFrame (cells x genes).')
    if max_cells is not None and int(max_cells) < 1:
        raise ValueError('max_cells must be None or a positive integer.')
    if int(predict_batch_size) < 1:
        raise ValueError('predict_batch_size must be a positive integer.')
    max_cells = None if max_cells is None else int(max_cells)
    predict_batch_size = int(predict_batch_size)

    if data_dir is None:
        raise ValueError('data_dir is required and must contain scFEA resource files.')
    data_dir = os.path.abspath(os.path.expanduser(data_dir))
    if not os.path.isdir(data_dir):
        raise FileNotFoundError('scFEA data_dir does not exist: %s' % data_dir)

    files = _SPECIES_FILES[species]
    moduleGene_file = files['moduleGene_file']
    cm_file = files['stoichiometry_matrix']
    cName_file = files['cName_file']

    # choose cpu or gpu automatically
    if device is None:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(device)

    # ------------------------------------------------------------------
    # read data (expression supplied directly as a cells x genes DataFrame)
    # ------------------------------------------------------------------
    if verbose:
        print("Starting load data...")
    geneExpr = expr.copy()
    geneExpr = geneExpr * 1.0
    if sc_imputation is True:
        import magic  # deferred import: only needed for imputation
        magic_operator = magic.MAGIC()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            geneExpr = magic_operator.fit_transform(geneExpr)
    if geneExpr.max().max() > 50:
        geneExpr = (geneExpr + 1).apply(np.log2)

    moduleGene = pd.read_csv(
        os.path.join(data_dir, moduleGene_file),
        sep=',',
        index_col=0)
    moduleLen = [moduleGene.iloc[i, :].notna().sum()
                 for i in range(moduleGene.shape[0])]
    moduleLen = np.array(moduleLen)

    # find existing gene
    module_gene_all = []
    for i in range(moduleGene.shape[0]):
        for j in range(moduleGene.shape[1]):
            if pd.isna(moduleGene.iloc[i, j]) is False:
                module_gene_all.append(moduleGene.iloc[i, j])
    module_gene_all = set(module_gene_all)
    data_gene_all = set(geneExpr.columns)
    gene_overlap = list(data_gene_all.intersection(module_gene_all))  # fix
    gene_overlap.sort()

    cmMat = pd.read_csv(
        os.path.join(data_dir, cm_file),
        sep=',',
        header=None)
    cmMat = cmMat.values
    cmMat = torch.FloatTensor(cmMat).to(device)

    cName = None
    if cName_file != 'noCompoundName':
        if verbose:
            print("Load compound name file, the balance output will have "
                  "compound name.")
        cName = pd.read_csv(
            os.path.join(data_dir, cName_file),
            sep=',',
            header=0)
        cName = cName.columns
    if verbose:
        print("Load data done.")

    if verbose:
        print("Starting process data...")
    # extract overlap gene
    geneExpr = geneExpr[gene_overlap]
    gene_names = list(geneExpr.columns)
    cell_names = geneExpr.index.astype(str).tolist()
    n_modules = moduleGene.shape[0]
    n_genes = len(gene_names)
    n_cells = len(cell_names)
    n_comps = cmMat.shape[0]

    # Pre-compute gene-name to column-index mapping
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    # Pre-compute per-module gene indices (and track empty modules)
    emptyNode = []
    module_gene_indices = []
    for i in range(n_modules):
        genes = moduleGene.iloc[i, :].values.astype(str)
        genes = [g for g in genes if g != 'nan']
        if not genes:
            emptyNode.append(i)
            module_gene_indices.append([])
            continue
        indices = [gene_to_idx[g] for g in genes if g in gene_to_idx]
        module_gene_indices.append(indices)

    # Convert expression to float32 numpy for fast tensor construction
    geneExpr_np = geneExpr.values.astype(np.float32)

    # Compute geneExprScale from the full dataset (needed for loss)
    geneExprScale = torch.from_numpy(
        geneExpr_np.sum(axis=1).astype(np.float32)
    ).float().to(device)
    stand = geneExprScale.mean()
    geneExprScale = geneExprScale / stand

    if verbose:
        print("Process data done.")

    # Determine training cells
    rng = np.random.RandomState(int(seed))
    all_indices = np.arange(n_cells, dtype=np.int64)
    if max_cells is not None and n_cells > max_cells:
        train_indices = rng.choice(all_indices, size=max_cells, replace=False)
        train_indices.sort()
        subset_msg = ("Training on %d randomly sampled cells "
                      "(max_cells=%d, total=%d)"
                      % (len(train_indices), max_cells, n_cells))
    else:
        train_indices = all_indices
        subset_msg = "Training on all %d cells" % n_cells
    if verbose:
        print(subset_msg)

    # Build training tensors only for the train subset
    X_train = _build_X_batch(geneExpr_np, module_gene_indices, n_genes,
                             n_modules, train_indices, device)
    module_scale_train = _build_module_scale(geneExpr_np, module_gene_indices,
                                             moduleLen, train_indices, device)
    geneExprScale_train = geneExprScale[train_indices]
    n_train = len(train_indices)

    # =====================================================================
    # NN
    torch.manual_seed(int(seed))
    net = FLUX(X_train, n_modules, f_in=n_genes, f_out=1).to(device)
    optimizer = torch.optim.Adam(net.parameters(), lr=LEARN_RATE)

    # Dataloader
    dataloader_params = {'batch_size': n_train,
                         'shuffle': False,
                         'num_workers': 0,
                         'pin_memory': False}

    dataSet = MyDataset(X_train, geneExprScale_train, module_scale_train)
    train_loader = torch.utils.data.DataLoader(dataset=dataSet,
                                               **dataloader_params)

    # =====================================================================
    if verbose:
        print("Starting train neural network...")
    start = time.time()
    # training
    loss_v = []
    loss_v1 = []
    loss_v2 = []
    loss_v3 = []
    loss_v4 = []
    net.train()
    epoch_iter = tqdm(range(n_epoch)) if verbose else range(n_epoch)
    for epoch in epoch_iter:
        loss, loss1, loss2, loss3, loss4 = 0, 0, 0, 0, 0

        for X_b, X_scale_b, m_scale_b in train_loader:

            X_batch = Variable(X_b.float().to(device))
            X_scale_batch = Variable(X_scale_b.float().to(device))
            m_scale_batch = Variable(m_scale_b.float().to(device))

            out_m_batch, out_c_batch = net(X_batch, n_modules, n_genes,
                                           n_comps, cmMat)
            (loss_batch, loss1_batch, loss2_batch, loss3_batch,
             loss4_batch) = myLoss(
                out_m_batch, out_c_batch,
                lamb1=LAMB_BA, lamb2=LAMB_NG, lamb3=LAMB_CELL, lamb4=LAMB_MOD,
                geneScale=X_scale_batch, moduleScale=m_scale_batch)

            optimizer.zero_grad()
            loss_batch.backward()
            optimizer.step()

            loss += loss_batch.cpu().data.numpy()
            loss1 += loss1_batch.cpu().data.numpy()
            loss2 += loss2_batch.cpu().data.numpy()
            loss3 += loss3_batch.cpu().data.numpy()
            loss4 += loss4_batch.cpu().data.numpy()

        loss_v.append(loss)
        loss_v1.append(loss1)
        loss_v2.append(loss2)
        loss_v3.append(loss3)
        loss_v4.append(loss4)

    # =====================================================================
    end = time.time()
    if verbose:
        print("Training time: ", end - start)

    # Prediction: batched inference over all cells
    if verbose:
        print("Starting prediction on %d cells (batch_size=%d)..."
              % (n_cells, predict_batch_size))

    fluxStatuTest = np.zeros((n_cells, n_modules), dtype=np.float32)
    balanceStatus = np.zeros((n_cells, n_comps), dtype=np.float32)
    net.eval()

    n_batches = int(np.ceil(n_cells / predict_batch_size))
    batch_iter = tqdm(range(n_batches)) if verbose else range(n_batches)
    with torch.no_grad():
        for bi in batch_iter:
            start_idx = bi * predict_batch_size
            end_idx = min(start_idx + predict_batch_size, n_cells)
            batch_indices = all_indices[start_idx:end_idx]

            X_batch = _build_X_batch(geneExpr_np, module_gene_indices,
                                     n_genes, n_modules, batch_indices, device)
            out_m, out_c = net(X_batch, n_modules, n_genes, n_comps, cmMat)
            fluxStatuTest[start_idx:end_idx, :] = out_m.detach().cpu().numpy()
            balanceStatus[start_idx:end_idx, :] = out_c.detach().cpu().numpy()

    # ------------------------------------------------------------------
    # assemble results as DataFrames (mirrors scFEA setF / setB outputs)
    # ------------------------------------------------------------------
    setF = pd.DataFrame(fluxStatuTest)
    setF.columns = moduleGene.index
    setF.index = geneExpr.index.tolist()

    setB = pd.DataFrame(balanceStatus)
    setB.rename(columns=lambda x: x + 1)
    setB.index = setF.index
    if cName_file != 'noCompoundName':
        setB.columns = cName

    if verbose:
        print("scFEA job finished.")

    predictions = {
        'flux': fluxStatuTest,
        'balance': balanceStatus,
        'cell_names': list(setF.index),
        'module_names': list(moduleGene.index),
        'compound_names': list(cName) if cName is not None else None,
        'loss': loss_v,
    }

    return {'flux': setF, 'balance': setB, 'predictions': predictions}
