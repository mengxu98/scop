"""Real Scanpy/PyNNDescent regression: python -m unittest discover -s tests/python."""
import importlib.metadata
import sys
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# Only the R thisutils logger is unavailable in standalone Python CI.
logger = types.ModuleType("log_message")
logger.log_message = lambda *args, **kwargs: None
sys.modules.setdefault("log_message", logger)
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "inst" / "python"))
from functions import PAGA, _check_paga_pynndescent_runtime


class PagaNeighborsTests(unittest.TestCase):
    def test_default_large_data_path(self):
        # Euclidean data below 8192 observations uses sklearn, hiding this bug.
        rng = np.random.default_rng(42)
        x = rng.normal(size=(8200, 10)).astype(np.float32)
        adata = sc.AnnData(sparse.csr_matrix(np.abs(x)))
        adata.obs["group"] = pd.Categorical(np.repeat(["A", "B", "C", "D"], 2050))
        adata.obsm["PCA"] = x
        adata.obsm["UMAP"] = x[:, :2].copy()
        with TemporaryDirectory() as directory:
            result = PAGA(
                adata=adata, group_by="group", linear_reduction="PCA",
                nonlinear_reduction="UMAP", n_pcs=10, n_neighbors=10,
                show_plot=False, dirpath=directory,
            )
        self.assertIs(result, adata)
        self.assertTrue(sparse.isspmatrix_csr(result.obsp["distances"]))
        graph = result.uns["paga"]["connectivities"].toarray()
        self.assertEqual(graph.shape, (4, 4))
        self.assertTrue(np.isfinite(graph).all())
        np.testing.assert_allclose(graph, graph.T)
        self.assertTrue(np.all((graph >= 0) & (graph <= 1)))
        self.assertGreater(np.count_nonzero(graph), 0)
        self.assertEqual(result.uns["paga"]["pos"].shape, (4, 2))
        self.assertIn("X_draw_graph_fr", result.obsm)

    def test_replaced_loaded_dependency_requests_restart(self):
        import pynndescent
        # Model pip replacing files while reticulate retains the loaded module.
        with patch.object(pynndescent, "__version__", "0.7.0"), patch(
            "importlib.metadata.version", return_value="0.6.0"
        ):
            with self.assertRaisesRegex(RuntimeError, "Restart R"):
                PAGA(adata=None)

    def test_matching_loaded_dependency_is_accepted(self):
        import pynndescent
        self.assertEqual(pynndescent.__version__, importlib.metadata.version("pynndescent"))
        _check_paga_pynndescent_runtime()


if __name__ == "__main__":
    unittest.main()
