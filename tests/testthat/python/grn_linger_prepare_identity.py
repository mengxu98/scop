import importlib.metadata
import importlib.util
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


runner_path = Path(sys.argv[1]).resolve()
spec = importlib.util.spec_from_file_location("linger_prepare_identity_test", runner_path)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
required = ("LingerGRN", "numpy", "pandas", "anndata", "scanpy", "torch", "shap")
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    (root / "runtime-source" / "LingerGRN").mkdir(parents=True)
    versions = {name: f"1.{index}" for index, name in enumerate(required)}
    for name, version in versions.items():
        metadata = root / f"{name}-{version}.dist-info"
        metadata.mkdir()
        (metadata / "METADATA").write_text(f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n")
    (root / "~nndata-0.9.2.dist-info").mkdir()
    unnamed = root / "unnamed-1.0.dist-info"
    unnamed.mkdir()
    (unnamed / "METADATA").write_text("Metadata-Version: 2.1\nVersion: 1.0\n")
    distributions = list(importlib.metadata.distributions(path=[str(root)]))
    request = dict(min_genes=20, min_cells=3, genome="hg38", activation="ReLU",
                   patch_missing_keys=True, seed=42, cores=1,
                   reference_selection={"manifest_sha256": "selected-spec-sha256"})
    with patch.object(runner.importlib.metadata, "distributions", return_value=distributions):
        valid = runner.prepare_identity(root, request, {"input": "same"}, {"reference": "same"}, {})
    assert valid["packages"] == versions
    assert valid["reference_spec_sha256"] == "selected-spec-sha256"
    with patch.object(runner.importlib.metadata, "distributions", return_value=distributions + [SimpleNamespace(metadata=None)]):
        malformed = runner.prepare_identity(root, request, {"input": "same"}, {"reference": "same"}, {})
    assert malformed == valid
    for name, version in versions.items():
        assert malformed["packages"][name] == version
print("LINGER_PREPARE_IDENTITY_METADATA_PASSED")
