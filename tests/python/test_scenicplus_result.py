import pickle
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
from scenicplus.scenicplus_mudata import ScenicPlusMuData

functions_file = Path(sys.argv.pop(1)).resolve()
export_dir = Path(sys.argv.pop(1)) if len(sys.argv) > 1 else None
sys.path.insert(0, str(functions_file.parent))
from functions import ExtractSCENICPlusTables


def write_fixture(directory, annotations=("direct", "extended")):
    directory.mkdir(parents=True, exist_ok=True)
    cells = ["c2", "c1", "c3"]
    frames = {}
    metadata = {}
    for annotation, tf, genes, regions, values in (
        ("direct", "TF1", ["G1", "G2"], ["r1", "r2"], [0.125, 0.5, 0.875]),
        ("extended", "TF2", ["G3", "G4"], ["r3", "r4"], [0.25, 0.625, 1.0]),
    ):
        if annotation not in annotations:
            continue
        signature = f"{tf}_{annotation}_+/+_(2g)"
        frames[annotation] = pd.DataFrame({signature: values}, index=cells)
        metadata[annotation] = pd.DataFrame({
            "TF": [tf, tf], "Gene": genes, "Region": regions,
            "Gene_signature_name": [signature, signature],
            "importance_TF2G": [0.75, 0.5], "importance_R2G": [0.625, 0.375],
        })
    auc = {}
    for annotation, frame in frames.items():
        if annotation == "extended":
            frame = frame.loc[["c3", "c2", "c1"]]
        auc[annotation] = md.MuData({
            "Gene_based": ad.AnnData(frame),
            "Region_based": ad.AnnData(frame.rename(columns=lambda x: x.replace("2g", "2r"))),
        })
    counts = md.MuData({
        "scRNA": ad.AnnData(pd.DataFrame(np.ones((3, 4)), index=cells, columns=["G1", "G2", "G3", "G4"])),
        "scATAC": ad.AnnData(pd.DataFrame(np.ones((3, 4)), index=cells, columns=["r1", "r2", "r3", "r4"])),
    })
    obj = ScenicPlusMuData(
        counts,
        e_regulon_auc_direct=auc.get("direct"),
        e_regulon_auc_extended=auc.get("extended"),
        e_regulon_metadata_direct=metadata.get("direct"),
        e_regulon_metadata_extended=metadata.get("extended"),
    )
    obj.write_h5mu(directory / "scplus.h5mu")
    combined = pd.concat(list(metadata.values()), ignore_index=True)
    reference = directory / "reference"
    reference.mkdir()
    for name, columns, rename in (
        ("tf_gene", ["TF", "Gene", "importance_TF2G"], {"Gene": "target", "importance_TF2G": "importance"}),
        ("region_gene", ["Region", "Gene", "importance_R2G"], {"Region": "region", "Gene": "gene", "importance_R2G": "score"}),
        ("triplets", ["TF", "Region", "Gene", "importance_R2G"], {"Region": "region", "Gene": "gene", "importance_R2G": "score"}),
    ):
        combined[columns].rename(columns=rename).to_csv(reference / f"{name}.tsv", sep="\t", index=False)
    tf_region = combined[["TF", "Region"]].rename(columns={"Region": "region"})
    tf_region["score"] = 1.0
    tf_region.to_csv(reference / "tf_region.tsv", sep="\t", index=False)
    eregulons = combined[["Gene_signature_name", "Gene"]].rename(columns={"Gene_signature_name": "regulon", "Gene": "target"})
    eregulons["target_count"] = 2
    eregulons.to_csv(reference / "eregulons.tsv", sep="\t", index=False)
    expected_auc = pd.concat(list(frames.values()), axis=1)
    expected_auc = expected_auc.loc[obj.mod[f"{annotations[0]}_gene_based_AUC"].obs_names]
    expected_auc.to_csv(reference / "auc.tsv", sep="\t")
    return obj


class ExtractSCENICPlusTests(unittest.TestCase):
    def test_reused_output_never_combines_different_objects(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            obj = write_fixture(directory, ("direct",))
            output = directory / "result"
            ExtractSCENICPlusTables(directory / "scplus.h5mu", output, verbose=False)
            previous = {p.name: p.read_bytes() for p in output.iterdir()}
            metadata = obj.uns["direct_e_regulon_metadata"].drop(columns="Region")
            with (directory / "incomplete.pkl").open("wb") as handle:
                pickle.dump({"uns": {"eRegulon_metadata": metadata, "eRegulon_AUC": {"Gene_based": obj.mod["direct_gene_based_AUC"].to_df()}}}, handle)
            with self.assertRaisesRegex(ValueError, "required SCENIC"):
                ExtractSCENICPlusTables(directory / "incomplete.pkl", output, verbose=False)
            self.assertEqual({p.name: p.read_bytes() for p in output.iterdir()}, previous)
            metadata = obj.uns["direct_e_regulon_metadata"].drop(columns="Gene_signature_name")
            with (directory / "without_signatures.pkl").open("wb") as handle:
                pickle.dump({"uns": {"eRegulon_metadata": metadata, "eRegulon_AUC": {"Gene_based": obj.mod["direct_gene_based_AUC"].to_df()}}}, handle)
            (output / "EREGULONS.CSV.gz").write_bytes(b"stale")
            (output / "keep.txt").write_text("unrelated file")
            ExtractSCENICPlusTables(directory / "without_signatures.pkl", output, verbose=False)
            self.assertFalse((output / "eregulons.tsv").exists())
            self.assertFalse((output / "EREGULONS.CSV.gz").exists())
            self.assertEqual((output / "keep.txt").read_text(), "unrelated file")

    def test_h5mu_preserves_dependency_import_errors(self):
        original_import = __import__

        def broken_import(name, *args, **kwargs):
            if name == "mudata":
                raise ImportError("controlled SciPy binary import failure")
            return original_import(name, *args, **kwargs)

        with tempfile.TemporaryDirectory() as tmp, patch("builtins.__import__", side_effect=broken_import):
            with self.assertRaisesRegex(ImportError, "controlled SciPy binary import failure"):
                ExtractSCENICPlusTables(Path(tmp) / "scplus.h5mu", Path(tmp) / "out", verbose=False)

    def test_publish_failure_restores_previous_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_fixture(directory, ("direct",))
            output = directory / "result"
            ExtractSCENICPlusTables(directory / "scplus.h5mu", output, verbose=False)
            previous = {p.name: p.read_bytes() for p in output.iterdir()}
            original_replace = Path.replace

            def fail_publish(path, target):
                if path.parent.name == "tables" and path.name == "tf_region.tsv":
                    raise OSError("controlled publication failure")
                return original_replace(path, target)

            with patch.object(Path, "replace", fail_publish):
                with self.assertRaisesRegex(OSError, "controlled publication failure"):
                    ExtractSCENICPlusTables(directory / "scplus.h5mu", output, verbose=False)
            self.assertEqual({p.name: p.read_bytes() for p in output.iterdir()}, previous)
            (output / "tf_gene.tsv").unlink()
            (output / "tf_gene.tsv").mkdir()
            with self.assertRaisesRegex(ValueError, "not a file"):
                ExtractSCENICPlusTables(directory / "scplus.h5mu", output, verbose=False)
            self.assertEqual((output / "auc.tsv").read_bytes(), previous["auc.tsv"])

    def test_runner_exports_without_duplicate_openmp_override(self):
        with tempfile.TemporaryDirectory(prefix="scenicplus result ") as tmp:
            directory = Path(tmp)
            obj = write_fixture(directory, ("direct",))
            with (directory / "legacy.pkl").open("wb") as handle:
                pickle.dump({"uns": {
                    "eRegulon_metadata": obj.uns["direct_e_regulon_metadata"],
                    "eRegulon_AUC": {"Gene_based": obj.mod["direct_gene_based_AUC"].to_df()},
                }}, handle)
            for name in ("scplus.h5mu", "legacy.pkl"):
                with self.subTest(name=name):
                    output = directory / f"{name} exported tables"
                    result = subprocess.run(
                        [sys.executable, str(functions_file.with_name("scenicplus_runner.py")),
                         str(directory / name), str(output), "--quiet"],
                        env={**os.environ, "KMP_DUPLICATE_LIB_OK": "FALSE"},
                        capture_output=True, text=True,
                    )
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    for reference in (directory / "reference").glob("*.tsv"):
                        expected = pd.read_csv(reference, sep="\t")
                        actual = pd.read_csv(output / reference.name, sep="\t")
                        pd.testing.assert_frame_equal(actual, expected)

    def test_h5mu_matches_official_direct_and_extended_results(self):
        for annotations in (("direct", "extended"), ("direct",), ("extended",)):
            with self.subTest(annotations=annotations), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                write_fixture(directory, annotations)
                ExtractSCENICPlusTables(directory / "scplus.h5mu", directory / "result", verbose=False)
                for reference in (directory / "reference").glob("*.tsv"):
                    index_col = 0 if reference.stem == "auc" else None
                    expected = pd.read_csv(reference, sep="\t", index_col=index_col)
                    actual = pd.read_csv(directory / "result" / reference.name, sep="\t", index_col=index_col)
                    pd.testing.assert_frame_equal(actual, expected)

    def test_legacy_pickle_keeps_gene_based_auc(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            obj = write_fixture(directory, ("direct",))
            with (directory / "legacy.pkl").open("wb") as handle:
                pickle.dump({"uns": {
                    "eRegulon_metadata": obj.uns["direct_e_regulon_metadata"],
                    "eRegulon_AUC": {"Gene_based": obj.mod["direct_gene_based_AUC"].to_df()},
                }}, handle)
            ExtractSCENICPlusTables(directory / "legacy.pkl", directory / "result", verbose=False)
            expected = pd.read_csv(directory / "reference" / "auc.tsv", sep="\t", index_col=0)
            actual = pd.read_csv(directory / "result" / "auc.tsv", sep="\t", index_col=0)
            pd.testing.assert_frame_equal(actual, expected)


if __name__ == "__main__":
    result = unittest.main(exit=False)
    if not result.result.wasSuccessful():
        sys.exit(1)
    if export_dir is not None:
        write_fixture(export_dir)
        ExtractSCENICPlusTables(export_dir / "scplus.h5mu", export_dir / "result", verbose=False)
