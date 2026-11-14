import importlib.metadata
import json
import os
import pickle
import signal
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from importlib.resources import files
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import yaml


RESOURCE_NAMES = (
    "cistopic", "gex", "region_sets", "ranking_db", "scores_db",
    "motif_annotation", "genome_annotation", "chromsizes",
)


class WorkflowCancelled(BaseException):
    def __init__(self, signum):
        self.signum = signum
        super().__init__("Official SCENIC+ cancelled by " + signal.Signals(signum).name)


@contextmanager
def execution_lock(output):
    try:
        import fcntl
    except ImportError as error:
        raise RuntimeError("The official SCENIC+ workflow requires Linux or macOS directory locking") from error
    path = output / ".scop-execution.lock"
    with path.open("a+b") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"Another official SCENIC+ process owns the execution directory: {output}") from error
        yield handle.fileno()


def atomic_json(path, value):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix="." + path.name + ".", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(value, handle, indent=2)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def run_workflow(command, output, environment, lock_fd):
    process = None
    cancelled = None
    previous = {}

    def cancel(signum, frame):
        nonlocal cancelled
        cancelled = signum

    try:
        for signum in (signal.SIGINT, signal.SIGTERM):
            previous[signum] = signal.signal(signum, cancel)
        with (output / "snakemake.stdout.log").open("w") as stdout, (output / "snakemake.stderr.log").open("w") as stderr:
            process = subprocess.Popen(command, cwd=output, env=environment, stdout=stdout, stderr=stderr,
                start_new_session=True, pass_fds=(lock_fd,))
            while True:
                if cancelled is not None:
                    raise WorkflowCancelled(cancelled)
                try:
                    returncode = process.wait(timeout=0.2)
                    if cancelled is not None:
                        raise WorkflowCancelled(cancelled)
                    return returncode
                except subprocess.TimeoutExpired:
                    continue
    finally:
        if process is not None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                pass
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def validate_request(request):
    resources = request["resources"]
    for name in RESOURCE_NAMES:
        path = Path(resources[name])
        if any(character.isspace() or character in "'\";$`\\|&<>(){}[]*?!" for character in str(path)):
            raise ValueError(f"The official SCENIC+ shell workflow requires resource paths without shell metacharacters: {name}")
    workflow = Path(resources["workflow"]) if "workflow" in resources else Path(str(files("scenicplus.snakemake")))
    if "workflow" not in resources and not all((workflow / name).is_file() for name in ("Snakefile", "config.yaml")):
        raise ValueError("The official workflow files are missing; set resources$workflow to the SCENIC+ snakemake source directory")
    return resources, workflow


def prepare_inputs(request, resources, output):
    cells = request["cells"]
    gex = ad.read_h5ad(resources["gex"])
    if not gex.obs_names.is_unique or not set(cells).issubset(gex.obs_names):
        raise ValueError("SCENIC+ GEX cell identifiers do not contain all requested cells")
    gex = gex[cells].copy()
    if request.get("expression"):
        values = scipy.io.mmread(request["expression"]).tocsr().T.tocsr()
        genes = Path(request["features"]).read_text().splitlines()
        if values.shape != (len(cells), len(genes)) or len(genes) != len(set(genes)):
            raise ValueError("SCENIC+ expression matrix dimensions or feature identifiers are invalid")
        metadata = gex.var.reindex(genes).copy()
        gex = ad.AnnData(X=values, obs=gex.obs.copy(), var=metadata)
    elif gex.raw is not None:
        gex = gex.raw.to_adata()
    if not np.isfinite(gex.X.data if hasattr(gex.X, "tocsr") else gex.X).all():
        raise ValueError("SCENIC+ expression contains nonfinite values")
    gex.raw = gex.copy()
    gex_path = output / "input_GEX.h5ad"
    gex.write_h5ad(gex_path)
    with open(resources["cistopic"], "rb") as handle:
        cistopic = pickle.load(handle)
    if len(cistopic.cell_names) != len(set(cistopic.cell_names)) or not set(cells).issubset(cistopic.cell_names):
        raise ValueError("SCENIC+ cisTopic cell identifiers do not contain all requested cells")
    if not hasattr(cistopic.selected_model, "cell_topic"):
        raise ValueError("SCENIC+ requires a cisTopic object with a selected fitted topic model")
    if set(cistopic.cell_names) != set(cells):
        index = [cistopic.cell_names.index(cell) for cell in cells]
        cistopic.fragment_matrix = cistopic.fragment_matrix[:, index]
        cistopic.binary_matrix = cistopic.binary_matrix[:, index]
        cistopic.cell_names = cells
        cistopic.cell_data = cistopic.cell_data.loc[cells].copy()
        cistopic.selected_model.cell_topic = cistopic.selected_model.cell_topic.loc[:, cells].copy()
        if "cell" in cistopic.projections:
            cistopic.projections["cell"] = {key: value.loc[cells].copy() for key, value in cistopic.projections["cell"].items()}
    cistopic_path = output / "input_cistopic.pkl"
    with cistopic_path.open("wb") as handle:
        pickle.dump(cistopic, handle, protocol=4)
    return str(cistopic_path), str(gex_path)


def configure(request, resources, workflow, output, prepared):
    config = yaml.safe_load((workflow / "config.yaml").read_text())
    config["input_data"].update(cisTopic_obj_fname=prepared[0], GEX_anndata_fname=prepared[1],
        region_set_folder=str(Path(resources["region_sets"]).resolve()), ctx_db_fname=str(Path(resources["ranking_db"]).resolve()),
        dem_db_fname=str(Path(resources["scores_db"]).resolve()), path_to_motif_annotations=str(Path(resources["motif_annotation"]).resolve()))
    config["params_general"].update(temp_dir=str(output / "tmp"), n_cpu=int(request["cores"]), seed=int(request["seed"]))
    config["params_data_preparation"].update(is_multiome=True, bc_transform_func='"lambda x: x"')
    config["output_data"]["genome_annotation"] = str(Path(resources["genome_annotation"]).resolve())
    config["output_data"]["chromsizes"] = str(Path(resources["chromsizes"]).resolve())
    args = dict(request["method_args"])
    args.pop("resume_dir", None)
    allowed = {"gsea_n_perm", "min_target_genes", "data_preparation", "motif_enrichment", "inference"}
    unknown = set(args) - allowed
    if unknown:
        raise ValueError("Unknown SCENIC+ method arguments: " + ", ".join(sorted(unknown)))
    for group in ("data_preparation", "motif_enrichment", "inference"):
        overrides = args.get(group, {})
        if not isinstance(overrides, dict):
            raise ValueError(f"SCENIC+ {group} must be a named configuration list")
        section = config["params_" + group]
        unknown = set(overrides) - set(section)
        if unknown:
            raise ValueError(f"Unknown SCENIC+ {group} arguments: " + ", ".join(sorted(unknown)))
        if group == "data_preparation" and set(overrides) & {"is_multiome", "bc_transform_func", "key_to_group_by"}:
            raise ValueError("SCENIC+ paired cell identity settings are controlled by the unified input")
        section.update(overrides)
    for name in ("gsea_n_perm", "min_target_genes"):
        if name in args:
            if name in args.get("inference", {}):
                raise ValueError(f"SCENIC+ argument {name} is specified twice")
            config["params_inference"][name] = args[name]
    path = output / "config.yaml"
    serialized = yaml.safe_dump(config, sort_keys=False)
    if not path.exists() or path.read_text() != serialized:
        path.write_text(serialized)
    return config, path


def export_result(output, cells):
    native = {}
    tables = {}
    identifiers = ("TF", "target", "region", "Region", "Gene",
                   "eRegulon_name", "Gene_signature_name", "Region_signature_name")
    for key, filename in dict(tf_gene="tf_to_gene_adj.tsv", region_gene="region_to_gene_adj.tsv",
        direct="eRegulon_direct.tsv", extended="eRegulons_extended.tsv").items():
        table = pd.read_csv(output / filename, sep="\t",
                            converters=dict.fromkeys(identifiers, str))
        table.to_csv(output / (key + ".tsv"), sep="\t", index=False, float_format="%.17g")
        native[key] = len(table)
        if key in ("direct", "extended"):
            tables[key] = table
    direct, extended = tables["direct"], tables["extended"]
    direct["annotation_type"] = "direct"
    extended["annotation_type"] = "extended"
    chains = pd.concat([direct, extended], ignore_index=True)
    chains.to_csv(output / "chains.tsv", sep="\t", index=False, float_format="%.17g")
    import scipy.sparse
    membership_parts = []
    for annotation in ("direct", "extended"):
        cistromes = ad.read_h5ad(output / ("cistromes_" + annotation + ".h5ad"))
        incidence = scipy.sparse.coo_matrix(cistromes.X)
        membership_parts.append(pd.DataFrame(dict(TF=np.asarray(cistromes.var_names)[incidence.col],
            Region=np.asarray(cistromes.obs_names)[incidence.row], annotation_type=annotation,
            membership=incidence.data.astype(np.float64))))
    memberships = pd.concat(membership_parts, ignore_index=True)
    memberships.to_csv(output / "tf_region.tsv", sep="\t", index=False)
    native["tf_region"] = len(memberships)
    native["chains"] = len(chains)
    gene_sets = {}
    for table in (direct, extended):
        if {"Gene_signature_name", "Gene"}.issubset(table):
            for signature, block in table.groupby("Gene_signature_name", sort=False):
                gene_sets[str(signature)] = block["Gene"].drop_duplicates().astype(str).tolist()
    (output / "gene_sets.json").write_text(json.dumps(gene_sets))
    import mudata
    obj = mudata.read_h5mu(output / "scplusmdata.h5mu")
    auc_files = []
    for name, modality in obj.mod.items():
        if "AUC" in name:
            table = modality.to_df()
            if set(table.index) != set(cells):
                raise ValueError(f"SCENIC+ AUC cells differ for {name}")
            table.loc[cells].to_csv(output / (name + ".tsv"), sep="\t", float_format="%.17g")
            auc_files.append(name + ".tsv")
    if not auc_files:
        raise ValueError("The completed SCENIC+ object has no AUC modalities")
    return native, auc_files


def main(request_path):
    request = json.loads(Path(request_path).read_text())
    resources, workflow = validate_request(request)
    output = Path(request["output_dir"]).resolve()
    if any(character.isspace() or character in "'\";$`\\|&<>(){}[]*?!" for character in str(output)):
        raise ValueError("The official SCENIC+ shell workflow requires a run directory without whitespace or shell metacharacters")
    output.mkdir(parents=True, exist_ok=True)
    with execution_lock(output) as lock_fd:
        run_locked(request, resources, workflow, output, lock_fd)


def run_locked(request, resources, workflow, output, lock_fd):
    manifest = output / "execution.json"
    resume = bool(request["method_args"].get("resume_dir"))
    if not resume and any(path.name != ".scop-execution.lock" for path in output.iterdir()):
        raise ValueError("Choose an empty SCENIC+ directory for a new run.")
    if manifest.exists() and json.loads(manifest.read_text())["status"] == "completed":
        raise ValueError("This SCENIC+ run is complete. Read scplusmdata.h5mu with RunSCENICPlus(scplus_object = ...) or select a new work_dir.")
    prepared = (str(output / "input_cistopic.pkl"), str(output / "input_GEX.h5ad"))
    marker = output / "prepare.completed"
    legacy_marker = output / "prepared_identity.json"
    if not (resume and (marker.is_file() or legacy_marker.is_file()) and
            all(Path(path).is_file() for path in prepared)):
        marker.unlink(missing_ok=True)
        legacy_marker.unlink(missing_ok=True)
        prepare_inputs(request, resources, output)
        marker.touch()
    config, config_path = configure(request, resources, workflow, output, prepared)
    (output / "tmp").mkdir(exist_ok=True)
    command = [str(Path(sys.executable).parent / "snakemake"), "--snakefile", str(workflow / "Snakefile"),
        "--configfile", str(config_path), "--cores", str(request["cores"]), "--rerun-incomplete", "--keep-going", "--latency-wait", "120", "--printshellcmds"]
    environment = os.environ.copy()
    environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
    environment["PYTHONHASHSEED"] = str(request["seed"])
    started = time.time()
    status = dict(status="running", implementation="SCENIC+::official_snakemake", command=command,
        input_mode="prepared_cistopic_paired_RNA", cells=request["cells"], parameters=config,
        package_versions={name: importlib.metadata.version(name) for name in ("scenicplus", "pycisTopic", "pycistarget", "snakemake")})
    atomic_json(manifest, status)
    try:
        returncode = run_workflow(command, output, environment, lock_fd)
        if returncode:
            raise RuntimeError(f"Official SCENIC+ exited {returncode}; see {output / 'snakemake.stderr.log'}")
        rows, auc_files = export_result(output, request["cells"])
        status.update(status="completed", elapsed_seconds=time.time() - started, rows=rows, auc_files=auc_files)
    except BaseException as error:
        status.update(status="cancelled" if isinstance(error, (WorkflowCancelled, KeyboardInterrupt)) else "failed",
            elapsed_seconds=time.time() - started, error=str(error))
        atomic_json(manifest, status)
        raise
    atomic_json(manifest, status)
    print(json.dumps(dict(status=status["status"], rows=status["rows"])), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
