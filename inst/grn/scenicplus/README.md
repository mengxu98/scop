# Run the official SCENIC+ workflow

The runner executes the official SCENIC+ Snakemake workflow from a fitted
cisTopic object, paired gene expression, candidate region sets and reference
databases. The fitted cisTopic model is an explicit input. Topic modelling and
reference database downloads are prepared separately.

Use `RunGRN(..., methods = "scenicplus", input_mode = "multiome",
return_data = TRUE)`. Supply these paths through `resources`.

| Resource | Input |
| --- | --- |
| `cistopic` | Pickled pycisTopic object with a selected fitted model |
| `gex` | Paired AnnData file containing the requested cell identifiers |
| `region_sets` | Directory of named collections containing BED region sets |
| `ranking_db` | cisTarget regions-versus-motifs rankings database |
| `scores_db` | DEM regions-versus-motifs scores database |
| `motif_annotation` | Motif-to-TF annotation table |
| `genome_annotation` | Offline SCENIC+ gene annotation table |
| `chromsizes` | Chromosome sizes table for that annotation |
| `workflow` | Optional official directory containing `Snakefile` and `config.yaml` |

The selected RNA assay and layer provide the expression values. The GEX file
provides paired cell metadata. Every requested cell must occur in both GEX and
cisTopic. The runner subsets the existing fitted topic posterior to the
requested cells and retains its original model identity. Use resources prepared
from the intended population when each population requires separately fitted
topics.

Some SCENIC+ distributions omit their workflow files. Set `workflow` explicitly
to the official source distribution's `src/scenicplus/snakemake` directory in
that case. The upstream workflow uses shell commands, so its resource paths
must have no whitespace or shell metacharacters.

Set `method_args$runtime$python` to the Python executable containing SCENIC+,
pycisTopic, pycistarget, Snakemake and pybiomart. The runner executes separate
processes. `cores` and `seed` configure the official workflow. Use
`method_args$inference`, `method_args$motif_enrichment` and
`method_args$data_preparation` for named overrides of the official configuration
sections. Unknown keys fail before execution. `gsea_n_perm` and
`min_target_genes` are convenient aliases for their inference keys.

The result retains the complete TF-gene and region-gene adjacency tables,
direct and extended eRegulon metadata, gene sets and every AUC modality. Each
adjacency layer uses its own correlation to describe its sign. TF-region
membership comes from the eRegulon cistrome and has an unknown regulatory sign.
The native `.h5mu` file, exact workflow configuration, execution receipt and
logs remain in the method run directory.

The official runner supports Linux and macOS. It holds an operating-system
lock on the resolved execution directory while preparing inputs, executing
Snakemake and exporting results. Concurrent calls with the same `resume_dir`
fail before changing that directory. Cancellation stops and waits for the
workflow subprocess group before releasing ownership. The lock file remains
as a stable file; the operating-system lock determines whether a run is active.
Request identities and execution receipts use atomic replacement.
Each output directory belongs to one scop workflow. A new directory must be
empty. Keep external native Snakemake workflows in separate directories.

To continue an interrupted run, set `method_args$resume_dir` to that run's
official output directory and keep its input content, resources, implementation
and scientific parameters unchanged. The runner validates those identities and
uses the official `--rerun-incomplete` workflow. Completed results are read from
the framework's verified cache.
After a hard interruption, automatic native-lock recovery requires the same
verified request identity, the exact lock files recorded after this workflow
started a job, and an inactive recorded process group. Recovery removes only
those recorded native lock files while retaining the scop directory lock.
It rechecks each file's content and inode before removal. Locks added by an
external workflow remain untouched and native locking still controls execution.
Do not run manual lock cleanup or replace files in this dedicated directory
while scop is recovering it.
Foreign, changed or unrecorded native locks require inspection and an official
manual unlock after their workflow has stopped. Native locking remains enabled.

Existing `RunSCENICPlus()` result import and C++ approximation keep their own
interfaces. This runner records `SCENIC+::official_snakemake` only after the
official workflow and result export complete.
