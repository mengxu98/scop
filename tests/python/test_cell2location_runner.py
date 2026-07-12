import importlib.util
import json
import sys
import types
from pathlib import Path

import pytest

ad = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")


RUNNER = Path(__file__).parents[2] / "inst" / "python" / "cell2location_runner.py"


def _load_runner():
    spec = importlib.util.spec_from_file_location("scop_cell2location_runner", RUNNER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _install_fake_backend(monkeypatch, seen):
    root = types.ModuleType("cell2location")
    root.__path__ = []
    models = types.ModuleType("cell2location.models")
    utils = types.ModuleType("cell2location.utils")
    utils.__path__ = []
    filtering = types.ModuleType("cell2location.utils.filtering")

    class RegressionModel:
        @classmethod
        def setup_anndata(cls, **kwargs):
            seen["reference_setup"] = kwargs

        def __init__(self, adata):
            self.adata = adata

        def train(self, **kwargs):
            seen["reference_train"] = kwargs

        def export_posterior(self, adata, sample_kwargs=None):
            seen["reference_posterior"] = sample_kwargs
            adata = adata.copy()
            adata.uns["mod"] = {"factor_names": ["Alpha", "Beta"]}
            adata.var["means_per_cluster_mu_fg_Alpha"] = [4.0, 1.0, 2.0]
            adata.var["means_per_cluster_mu_fg_Beta"] = [1.0, 4.0, 2.0]
            return adata

        def save(self, path, overwrite=False):
            Path(path).mkdir(parents=True, exist_ok=True)
            (Path(path) / "model.pt").write_text("reference", encoding="utf-8")

    class Cell2location:
        @classmethod
        def setup_anndata(cls, **kwargs):
            seen["spatial_setup"] = kwargs

        def __init__(self, adata, **kwargs):
            self.adata = adata
            seen["spatial_init"] = kwargs

        def train(self, **kwargs):
            seen["spatial_train"] = kwargs

        def export_posterior(self, adata, **kwargs):
            seen["spatial_posterior"] = kwargs
            adata = adata.copy()
            adata.uns["mod"] = {"factor_names": ["Alpha", "Beta"]}
            adata.obsm["q05_cell_abundance_w_sf"] = pd.DataFrame(
                [[8.0, 2.0], [2.0, 6.0]],
                index=adata.obs_names,
                columns=["Alpha", "Beta"],
            )
            return adata

        def save(self, path, overwrite=False):
            Path(path).mkdir(parents=True, exist_ok=True)
            (Path(path) / "model.pt").write_text("spatial", encoding="utf-8")

    models.RegressionModel = RegressionModel
    models.Cell2location = Cell2location
    filtering.filter_genes = lambda adata, **kwargs: np.ones(adata.n_vars, dtype=bool)
    monkeypatch.setitem(sys.modules, "cell2location", root)
    monkeypatch.setitem(sys.modules, "cell2location.models", models)
    monkeypatch.setitem(sys.modules, "cell2location.utils", utils)
    monkeypatch.setitem(sys.modules, "cell2location.utils.filtering", filtering)


def _write_inputs(tmp_path):
    reference = ad.AnnData(
        X=np.array([[4, 1, 2], [1, 4, 2], [3, 1, 2], [1, 3, 2]], dtype=np.float32),
        obs=pd.DataFrame(
            {"celltype": ["Alpha", "Beta", "Alpha", "Beta"], "sample": ["r1", "r1", "r2", "r2"]},
            index=["Cell1", "Cell2", "Cell3", "Cell4"],
        ),
        var=pd.DataFrame(index=["Gene1", "Gene2", "Gene3"]),
    )
    spatial = ad.AnnData(
        X=np.array([[5, 2, 1], [1, 4, 2]], dtype=np.float32),
        obs=pd.DataFrame({"sample": ["s1", "s1"]}, index=["Spot1", "Spot2"]),
        var=pd.DataFrame(index=["Gene1", "Gene2", "Gene3"]),
    )
    reference_path = tmp_path / "reference.h5ad"
    spatial_path = tmp_path / "spatial.h5ad"
    reference.write_h5ad(reference_path)
    spatial.write_h5ad(spatial_path)
    return reference_path, spatial_path


def test_two_stage_runner_and_resume(tmp_path, monkeypatch):
    seen = {}
    _install_fake_backend(monkeypatch, seen)
    runner = _load_runner()
    monkeypatch.setattr(runner, "_versions", lambda: {"cell2location": "0.1.5", "scvi-tools": "1.3.3"})
    reference_path, spatial_path = _write_inputs(tmp_path)
    result_dir = tmp_path / "result with spaces"
    config = {
        "result_dir": str(result_dir),
        "spatial_path": str(spatial_path),
        "reference_path": str(reference_path),
        "signatures_path": None,
        "reference_label": "celltype",
        "spatial_batch": "sample",
        "reference_batch": "sample",
        "reference_covariates": None,
        "N_cells_per_location": 30,
        "detection_alpha": 20,
        "gene_filter_params": {"cell_count_cutoff": 1},
        "reference_train_params": {"max_epochs": 1, "accelerator": "cpu", "devices": "auto"},
        "spatial_train_params": {"max_epochs": 2, "accelerator": "cpu", "device": "auto"},
        "reference_posterior_params": {"num_samples": 2, "batch_size": 2},
        "spatial_posterior_params": {"batch_size": 2},
        "resume": True,
        "overwrite": False,
    }

    runner.run(config.copy())
    abundance = pd.read_csv(result_dir / "tables" / "abundance_q05.csv", index_col=0)
    proportions = pd.read_csv(result_dir / "tables" / "proportions.csv", index_col=0)
    assert list(abundance.columns) == ["Alpha", "Beta"]
    np.testing.assert_allclose(proportions.sum(axis=1), 1.0)
    assert seen["reference_train"]["accelerator"] == "cpu"
    assert seen["reference_train"]["device"] == "auto"
    assert "devices" not in seen["reference_train"]
    assert seen["spatial_train"]["max_epochs"] == 2
    assert seen["spatial_posterior"]["use_quantiles"] is True
    assert (result_dir / "reference" / "model" / "model.pt").exists()
    assert (result_dir / "spatial" / "model" / "model.pt").exists()

    seen.clear()
    runner.run(config.copy())
    assert seen == {}
    manifest = json.loads((result_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["reused_reference"] is True
    assert manifest["reused_spatial"] is True


def test_runner_rejects_mismatched_resume(tmp_path, monkeypatch):
    seen = {}
    _install_fake_backend(monkeypatch, seen)
    runner = _load_runner()
    monkeypatch.setattr(runner, "_versions", lambda: {})
    _, spatial_path = _write_inputs(tmp_path)
    signatures = pd.DataFrame(
        [[4.0, 1.0], [1.0, 4.0], [2.0, 2.0]],
        index=["Gene1", "Gene2", "Gene3"],
        columns=["Alpha", "Beta"],
    )
    signatures_path = tmp_path / "signatures.csv"
    signatures.to_csv(signatures_path)
    result_dir = tmp_path / "result"
    config = {
        "result_dir": str(result_dir),
        "spatial_path": str(spatial_path),
        "reference_path": None,
        "signatures_path": str(signatures_path),
        "reference_label": None,
        "spatial_batch": "sample",
        "reference_batch": None,
        "reference_covariates": None,
        "N_cells_per_location": 30,
        "detection_alpha": 20,
        "gene_filter_params": {},
        "reference_train_params": {},
        "spatial_train_params": {"max_epochs": 1},
        "reference_posterior_params": {},
        "spatial_posterior_params": {},
        "resume": True,
        "overwrite": False,
    }
    runner.run(config.copy())
    config["detection_alpha"] = 200
    with pytest.raises(RuntimeError, match="do not match"):
        runner.run(config.copy())
