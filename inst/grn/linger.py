from __future__ import annotations

import csv
import gzip
import importlib.metadata
import importlib.util
import json
import os
import random
import shutil
import sys
import tempfile
import time
from pathlib import Path


PREPARE_FILES = ("data/adata_RNA.h5ad", "data/adata_ATAC.h5ad", "data/TG_pseudobulk.tsv", "data/RE_pseudobulk.tsv", "data/Peaks.txt", "retained_cells.tsv", "input-summary.json")


def write_prepare_checkpoint(root):
    candidate = root / ".prepare.completed.json"
    candidate.write_text('{"status": "completed"}\n')
    candidate.replace(root / "prepare.completed.json")


def restore_prepare_checkpoint(source, root):
    source = Path(source).resolve()
    if json.loads((source / "prepare.completed.json").read_text()).get("status") != "completed":
        raise ValueError("LINGER resume_dir must contain a completed prepare checkpoint.")
    for name in PREPARE_FILES:
        destination = root / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / name, destination)
    write_prepare_checkpoint(root)
    return {"stage": "prepare", "source": str(source)}


def canonical_region(region):
    if ":" in region:
        chromosome, interval = region.rsplit(":", 1)
        start, end = interval.split("-", 1)
    else:
        chromosome, start, end = region.rsplit("-", 2)
    start, end = int(start), int(end)
    if start < 0 or end <= start:
        raise ValueError(f"LINGER requires valid primary-chromosome regions: {region}")
    return f"{chromosome}:{start}-{end}"


def prepare_dependency(root, patch_missing_keys):
    binary = str(Path(sys.executable).parent)
    os.environ["PATH"] = os.pathsep.join([binary] + [path for path in os.environ.get("PATH", "").split(os.pathsep) if path != binary])
    if shutil.which("bedtools") is None:
        raise ValueError("LINGER requires bedtools in the selected Python environment or PATH.")
    temporary = root / "tmp"
    temporary.mkdir()
    os.environ["TMPDIR"] = str(temporary)
    tempfile.tempdir = str(temporary)
    import pybedtools
    pybedtools.set_tempdir(str(temporary))
    version = importlib.metadata.version("LingerGRN")
    if version != "1.110":
        raise ValueError(f"LINGER requires LingerGRN 1.110, found {version}")
    source = Path(importlib.util.find_spec("LingerGRN.LL_net").origin).parent
    network_source = (source / "LL_net.py").read_text()
    before = "        for k in range(len(TFidxtemp)):\n            TG_1.append(geneName[ii])\n            TF_1.append(TFName_temp[k])\n            score_1.append(zscored_arr[k])"
    after = "            for k in range(len(TFidxtemp)):\n                TG_1.append(geneName[ii])\n                TF_1.append(TFName_temp[k])\n                score_1.append(zscored_arr[k])"
    function = network_source.split("def trans_shap(", 1)[1].split("\ndef ", 1)[0]
    needed = before in function
    if needed and not patch_missing_keys:
        raise ValueError("LingerGRN 1.110 requires the trans_shap missing-key patch; enable patch_missing_keys.")
    if not needed and after not in function:
        raise ValueError("Unrecognized LingerGRN trans_shap implementation; expected version 1.110 and its documented patch.")
    small_before = "    resultlist=[0 for i in range(times+1)]\n    for ii in tqdm(range(times), desc=chr):"
    small_after = "    resultlist=[0 for i in range(times+1)]\n    ii=-1\n    for ii in tqdm(range(times), desc=chr):"
    small_function = network_source.split("def TF_RE_LINGER_chr(", 1)[1].split("\ndef ", 1)[0]
    small_needed = small_before in small_function
    if not small_needed and small_after not in small_function:
        raise ValueError("Unrecognized LingerGRN TF_RE_LINGER_chr batch implementation.")
    destination = root / "runtime-source" / "LingerGRN"
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    license_source = Path(importlib.metadata.distribution("LingerGRN").locate_file("lingergrn-1.110.dist-info/LICENSE"))
    if license_source.is_file():
        shutil.copyfile(license_source, destination / "LICENSE")
    patched_source = network_source.replace(before, after, 1) if needed else network_source
    if small_needed:
        patched_source = patched_source.replace(small_before, small_after, 1)
    (destination / "LL_net.py").write_text(patched_source)
    preprocessing_source = (destination / "preprocess.py").read_text()
    preprocessing = preprocessing_source.replace("n_jobs = 8", "n_jobs = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))")
    (destination / "preprocess.py").write_text(preprocessing)
    sys.path.insert(0, str(destination.parent))
    sys.modules.pop("LingerGRN", None)
    importlib.invalidate_caches()
    return {
        "package_version": version, "missing_key_patch": "trans_shap-missing-key-v1",
        "patch_applied": needed, "worker_limit": int(os.environ["SLURM_CPUS_PER_TASK"]),
        "small_batch_patch": "tf-re-small-batch-v1", "small_batch_patch_applied": small_needed,
    }


def prepare(root, request):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy.io import mmread
    from LingerGRN.pseudo_bulk import pseudo_bulk

    inputs = Path(request["input_dir"])
    metadata = pd.read_csv(inputs / "cells.tsv", sep="\t", dtype=str, keep_default_na=False)
    if list(metadata.columns) != ["cell_id", "sample", "label"] or metadata.empty:
        raise ValueError("Cell metadata must contain cell_id, sample and label.")
    if metadata.isin([""]).any().any() or metadata.cell_id.duplicated().any():
        raise ValueError("LINGER cell metadata must have unique cells and complete sample/label values.")
    metadata = metadata.set_index("cell_id", drop=False)
    metadata["barcode"] = metadata.index
    genes = (inputs / "genes.txt").read_text().splitlines()
    regions = [canonical_region(region) for region in (inputs / "regions.txt").read_text().splitlines()]
    if len(set(genes)) != len(genes) or len(set(regions)) != len(regions):
        raise ValueError("LINGER feature identifiers must be unique.")
    rna = mmread(inputs / "rna.mtx").tocsc()
    atac = mmread(inputs / "atac.mtx").tocsc()
    if rna.shape != (len(genes), len(metadata)) or atac.shape != (len(regions), len(metadata)):
        raise ValueError("LINGER matrix dimensions do not match cell and feature identifiers.")
    rna_data = ad.AnnData(rna.T.tocsr(), obs=metadata.copy(), var=pd.DataFrame(index=genes))
    atac_data = ad.AnnData(atac.T.tocsr(), obs=metadata.copy(), var=pd.DataFrame(index=regions))
    rna_data.var["gene_ids"] = genes
    atac_data.var["gene_ids"] = regions
    for data in (rna_data, atac_data):
        sc.pp.filter_cells(data, min_genes=request["min_genes"])
        sc.pp.filter_genes(data, min_cells=request["min_cells"])
    common = [cell for cell in metadata.index if cell in rna_data.obs_names and cell in atac_data.obs_names]
    rna_data, atac_data = rna_data[common].copy(), atac_data[common].copy()
    if len(common) < 50:
        raise ValueError("LINGER filtering retained fewer than 50 paired cells.")
    data_dir = root / "data"
    data_dir.mkdir()
    rna_data.write_h5ad(data_dir / "adata_RNA.h5ad", compression="gzip")
    atac_data.write_h5ad(data_dir / "adata_ATAC.h5ad", compression="gzip")
    random.seed(request["seed"])
    np.random.seed(request["seed"])
    samples = sorted(rna_data.obs["sample"].unique())
    label_counts = rna_data.obs.groupby(["sample", "label"], observed=True).size()
    estimated = sum((1 if len(samples) ** 2 > 100 else int(np.floor(np.sqrt(count))) + 1) for count in label_counts if count >= 10)
    if estimated < 50:
        raise ValueError(f"LINGER sample/label groups can produce only {estimated} native pseudobulks; training requires at least 50.")
    tg_parts, re_parts = [], []
    for sample in samples:
        selected = rna_data.obs["sample"] == sample
        if selected.sum() < 20:
            raise ValueError(f"LINGER sample {sample} has fewer than 20 paired cells for native PCA/neighbors.")
        tg, re = pseudo_bulk(rna_data[selected].copy(), atac_data[selected].copy(), len(samples) ** 2 > 100)
        tg_parts.append(tg)
        re_parts.append(re)
    tg = pd.concat(tg_parts, axis=1).fillna(0)
    re = pd.concat(re_parts, axis=1).fillna(0).clip(upper=100)
    if tg.columns.duplicated().any() or list(tg.columns) != list(re.columns):
        raise ValueError("LINGER pseudobulk identities are duplicated or modalities are unaligned.")
    if tg.shape[1] < 50:
        raise ValueError(f"LINGER produced {tg.shape[1]} pseudobulks; native training requires at least 50 for its fixed SHAP background.")
    pd.Series(re.index).to_csv(data_dir / "Peaks.txt", sep="\t", header=False, index=False)
    tg.to_csv(data_dir / "TG_pseudobulk.tsv", sep="\t")
    re.to_csv(data_dir / "RE_pseudobulk.tsv", sep="\t")
    metadata.loc[common].to_csv(root / "retained_cells.tsv", sep="\t", index=False)
    summary = {"source_cells": len(metadata), "retained_cells": len(common), "samples": len(samples), "rna_features": rna_data.n_vars, "atac_features": atac_data.n_vars, "pseudobulks": tg.shape[1]}
    (root / "input-summary.json").write_text(json.dumps(summary, indent=2))


def preprocess(root, request):
    import pandas as pd
    from LingerGRN.preprocess import preprocess as native_preprocess

    data_dir, output = root / "data", root / "output"
    output.mkdir()
    tg = pd.read_csv(data_dir / "TG_pseudobulk.tsv", sep="\t", index_col=0)
    re = pd.read_csv(data_dir / "RE_pseudobulk.tsv", sep="\t", index_col=0)
    if request.get("regulators") or request.get("targets"):
        candidates = set(request.get("regulators", [])) | set(request.get("targets", []))
        tg = tg.loc[tg.index.isin(candidates)]
        if tg.empty:
            raise ValueError("No requested LINGER candidate genes remain in the native pseudobulk matrix.")
    native_preprocess(tg, re, request["bulk_reference"] + "/", request["genome"], "LINGER", str(output) + "/")
    required = ["Exp.txt", "Symbol.txt", "TFexp.txt", "TFName.txt", "Openness.txt", "TF_binding.txt", "index.txt"]
    missing = [name for name in required if not (output / name).is_file()]
    if missing:
        raise ValueError("LINGER preprocessing did not produce: " + ", ".join(missing))


def train(root, request):
    from LingerGRN import LINGER_tr

    LINGER_tr.training(request["bulk_reference"] + "/", "LINGER", str(root / "output") + "/", request["activation"], "Human")
    validate_training(root)


def validate_training(root):
    import pandas as pd
    import torch

    output = root / "output"
    genes = pd.read_csv(output / "data_merge.txt", sep="\t", index_col=0)
    chromosomes = sorted(pd.read_csv(output / "Region.bed", sep="\t", header=None)[0].unique())
    rows = []
    for chromosome in chromosomes:
        expected = set(map(int, genes.index[genes["chr"] == chromosome]))
        net = set(map(int, torch.load(output / f"net_{chromosome}.pt", map_location="cpu", weights_only=False)))
        shap = set(map(int, torch.load(output / f"shap_{chromosome}.pt", map_location="cpu", weights_only=False)))
        result = pd.read_csv(output / f"result_{chromosome}.txt", sep="\t", index_col=0)
        failed = set(result.index[pd.to_numeric(result.iloc[:, 0], errors="coerce") == -100])
        skipped = expected - net
        skipped_symbols = set(genes.loc[list(skipped), "Symbol"].astype(str)) if skipped else set()
        if net != shap or net - expected or skipped_symbols != failed:
            raise ValueError(f"Inconsistent LINGER trained models, SHAP and native skip records for {chromosome}.")
        rows.append(dict(chromosome=chromosome, expected_genes=len(expected), trained_genes=len(net), method_skipped_genes=len(skipped), method_skipped_symbols=",".join(sorted(skipped_symbols))))
    pd.DataFrame(rows).to_csv(root / "training-coverage.tsv", sep="\t", index=False)
    if not sum(row["trained_genes"] for row in rows):
        raise ValueError("LINGER did not train any target model.")


def networks(root, request):
    import anndata as ad
    from LingerGRN import LL_net

    output = str(root / "output") + "/"
    rna = ad.read_h5ad(root / "data/adata_RNA.h5ad")
    atac = ad.read_h5ad(root / "data/adata_ATAC.h5ad")
    reference = request["bulk_reference"] + "/"
    LL_net.TF_RE_binding(reference, rna, atac, request["genome"], "LINGER", output)
    LL_net.cis_reg(reference, rna, atac, request["genome"], "LINGER", output)
    LL_net.trans_reg(reference, "LINGER", output, request["genome"])


def write_dense_edges(source, destination, row_name, column_name, sources=None, targets=None):
    import numpy as np
    import pandas as pd

    count = 0
    with gzip.open(destination, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_name, column_name, "score"])
        for chunk in pd.read_csv(source, sep="\t", index_col=0, chunksize=500, dtype=str):
            for (row, column), score in chunk.stack(dropna=True).items():
                if sources is not None and str(column) not in sources:
                    continue
                if targets is not None and str(row) not in targets:
                    continue
                value = float(score)
                if np.isfinite(value) and value != 0:
                    writer.writerow([row, column, f"{value:.17g}"])
                    count += 1
    return count


def export(root, request):
    import numpy as np
    import pandas as pd

    output, network_dir = root / "output", root / "networks"
    network_dir.mkdir()
    regulators = set(request["regulators"]) if request.get("regulators") else None
    targets = set(request["targets"]) if request.get("targets") else None
    counts = {}
    counts["tf_region"] = write_dense_edges(output / "cell_population_TF_RE_binding.txt", network_dir / "tf_region.tsv.gz", "region", "regulator", sources=regulators)
    counts["tf_gene"] = write_dense_edges(output / "cell_population_trans_regulatory.txt", network_dir / "tf_gene.tsv.gz", "target", "regulator", sources=regulators, targets=targets)
    count = 0
    with gzip.open(network_dir / "region_gene.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["region", "target", "score"])
        for chunk in pd.read_csv(output / "cell_population_cis_regulatory.txt", sep="\t", header=None, names=["region", "target", "score"], dtype={"region": str, "target": str}, chunksize=500000):
            chunk["score"] = pd.to_numeric(chunk["score"], errors="coerce")
            chunk = chunk[np.isfinite(chunk.score) & (chunk.score != 0)]
            if targets is not None:
                chunk = chunk[chunk.target.isin(targets)]
            for row in chunk.itertuples(index=False):
                writer.writerow([row.region, row.target, f"{float(row.score):.17g}"])
                count += 1
    counts["region_gene"] = count
    pd.DataFrame({"endpoint": list(counts), "relations": list(counts.values())}).to_csv(network_dir / "network_summary.tsv", sep="\t", index=False)
    return counts


def run(request):
    root = Path(request["output_dir"]).resolve()
    if root.exists() and any(root.iterdir()):
        raise ValueError(f"LINGER output directory must be empty: {root}")
    root.mkdir(parents=True, exist_ok=True)
    inputs = Path(request["input_dir"]).resolve()
    request["input_dir"] = str(inputs)
    request["bulk_reference"] = str(Path(request["bulk_reference"]).resolve())
    chromosomes = sorted({canonical_region(region).split(":")[0] for region in (inputs / "regions.txt").read_text().splitlines()})
    selection = request["reference_selection"]
    if selection["genome"] != request["genome"] or selection["chromosomes"] != chromosomes:
        raise ValueError("LINGER reference selection does not match the requested genome and peak chromosomes.")
    reference_files = request["expected_reference_files"]
    if not reference_files or any(Path(name).name != name or name in (".", "..") for name in reference_files):
        raise ValueError("LINGER requires selected bulk reference file names.")
    resume_dir = request.get("resume_dir")
    checkpoint = request.get("checkpoint_prepare", True) or bool(resume_dir)
    dependency = prepare_dependency(root, request["patch_missing_keys"])
    resume = restore_prepare_checkpoint(resume_dir, root) if resume_dir else None
    before = Path.cwd()
    timing = []
    os.chdir(root)
    try:
        for step in (prepare, preprocess, train, networks, export):
            phase = step.__name__
            if step is prepare and resume is not None:
                timing.append({"phase": phase, "seconds": 0.0, "reused": True})
                print("phase_reused\tprepare", flush=True)
                continue
            print(f"phase_started\t{phase}", flush=True)
            started = time.perf_counter()
            result = step(root, request)
            if step is prepare and checkpoint:
                write_prepare_checkpoint(root)
            timing.append({"phase": phase, "seconds": time.perf_counter() - started})
            print(f"phase_finished\t{phase}", flush=True)
    finally:
        os.chdir(before)
    receipt = {
        "status": "completed", "implementation": "LingerGRN::LINGER[1.110;trans_shap-missing-key-v1;tf-re-small-batch-v1]",
        "dependency": dependency, "parameters": request, "phase_timings": timing,
        "reference_files": reference_files, "reference_selection": selection,
        "networks": result,
        "resume": resume,
        "versions": {name: importlib.metadata.version(name) for name in ("LingerGRN", "numpy", "pandas", "anndata", "scanpy", "torch", "shap")},
    }
    candidate = root / ".receipt.json"
    candidate.write_text(json.dumps(receipt, indent=2))
    candidate.replace(root / "receipt.json")
    return receipt


if __name__ == "__main__":
    request = json.loads(Path(sys.argv[1]).read_text())
    os.environ["SLURM_CPUS_PER_TASK"] = str(request["cores"])
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    run(request)
