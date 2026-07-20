# scop 0.9.0

* **feat**:
  * Added `RunBenchmark()` for isolated, failure-tolerant comparison of stable spatial-domain clustering methods against a gold standard, with ARI/NMI/purity, runtime, and sampled process-tree peak memory. Explicit `GiottoCluster` and `methods = "all"` selection add the supported legacy baseline without mixing it into the stable default, and result tables/plots label its tier. `BenchmarkPlot()` now provides publication-oriented quality, efficiency, overview, and direction-aware heatmap views for `scop_benchmark` objects.
  * Fixed Giotto runtime checks to use the current `drieslab/Giotto` repository registered by `SpatialBackendStatus()`, eliminating a false available-in-registry but unavailable-at-execution state.
  * `RunSpatialNeighborhood()` now defaults to truthful native observed summaries; `method = "spicyR"` requires `split.by` and records spicyR provenance only after the backend runs. The standard spatial workflow now records per-stage requested/completed/skipped/failed state and reports partial completion when a requested stage is skipped.
  * Registered `RunCell2location()` as a stable Python spatial deconvolution producer with schema-v1 results, read-only environment diagnostics in `SpatialBackendStatus()`, and bidirectional producer registry contracts.
  * Spatial analysis, plotting, and framework conversion now share strict image selection: objects with multiple spatial images must provide `image`, preventing silent first-slice truncation. Fifteen distance-sensitive producers now default to raw acquisition coordinates; their explicit `coordinate_space = "legacy_display"` path remains available for reproducing older display-scaled analyses.
  * Added `RunSpatialCellChat()` as a runtime-optional SpatialCellChat v3 producer for cell-, spot-, and composition-level communication. It enforces raw-to-micron coordinate contracts, isolates samples before distance calculations, stores truthful schema-v1 provenance, keeps native objects optional, and publishes group-level results under the distinct `SpatialCellChat` CCC method. `GetCCCObject()` and `SpatialCellChatPlot()` provide native-object access and stored-result spatial views without rerunning the backend.
  * Added spatial contract schema v1 with `SpatialCoordinates()`, `GetSpatialResult()`, and `GetSpatialGraph()`. Small-method results now retain their existing top-level payload fields while adding shared source, provenance, parameter, summary, and schema metadata; legacy framework results are normalized only in read-only views.
  * Added a shared spatial method registry with `ListSpatialMethods()`, read-only `SpatialBackendStatus()`, and `SpatialResultInfo()` so users can discover the complete spatial API, inspect optional backend compatibility without installation side effects, and summarize legacy or current results stored in `@tools`.
  * Added a strict raw-coordinate contract and sparse spatial graph core. `RunSpatialNetwork()` preserves raw distances and deterministic graph slots, and distance-sensitive producers expose raw coordinates as the analysis default while retaining an explicit legacy display path.
  * Replaced dense full-distance matrices in SmoothClust, spatial-variable-feature, and spatial-neighborhood paths with `BiocNeighbors` queries, including deterministic handling for duplicate coordinates and valid zero-edge radius graphs.
  * Added `MistyRPlot()` and `StatialKontextualPlot()` to close the visualization loop for stored high-level spatial context results without wrapping or rerunning their native frameworks.
  * Added `RunCell2location()` for official Python cell2location reference-signature learning and spatial abundance mapping through a persisted subprocess workflow. The wrapper records q05 absolute abundance and normalized proportions, supports manifest-validated resume and signature reuse, integrates with `standard_scop()`, and adds `Cell2locationPlot()` for abundance, proportion, dominant-type, and pie maps. `PrepareEnv(modules = "cell2location")` installs cell2location 0.1.5 with scvi-tools 1.3.3 in the shared SCOP environment.
  * `FeatureStatPlot()` and `ExpressionStatPlot()`: Added `auto_comparison` to automatically compare the group with the highest median statistic against all other groups, with explicit `ref_group` and `comparisons` still taking precedence.
  * `GSVAPlot()` now supports `mode = "diff"` for true two-group pathway activity tests on sample-aggregated GSVA scores, while keeping the original score plotting behavior as the default. Non-heatmap score-mode p-value columns are documented as score-derived plotting placeholders, not statistical significance.
  * Added `ClusterTreePlot()` for SCOP-styled visualization of Seurat multi-resolution clustering trees, including automatic `*_snn_res.*` metadata detection, prefix/resolution filtering, edge contribution statistics, and marker-expression overlays.
  * Added `RunCNV()` as a unified CNV workflow for Seurat objects, with runtime-optional `copykat`, `infercnv`, `SCEVAN`, and `fastCNV` backends, standardized result storage in `srt@tools[["CNV"]]`, metadata writeback, and `CNVPlot()` heatmap, embedding, spatial, composition-bar, and tree visualizations.
  * Added `RunESTIMATE()` for tumor microenvironment scoring from Seurat, `SummarizedExperiment`, or matrix inputs, including stromal, immune, ESTIMATE, and tumor-purity scores, Seurat metadata writeback, and `EstimateScorePlot()` / `EstimateGenePlot()` visualization helpers.
  * `RunMetabolism()`: Gene sets are now built via `PrepareDB()` by default (`use_preparedb = TRUE`) for species-aware gene mapping through BioMart and KEGG/Reactome databases. The `species` parameter now automatically converts human gene symbols to the target species. scMetabolism-curated pathway lists are cross-referenced with PrepareDB TERM2GENE so mouse data receives mouse gene symbols directly. The previous GMT-only path is still available with `use_preparedb = FALSE`.
  * `RunMetabolism()`: `convert_species` now defaults to `TRUE`, enabling automatic `GeneConvert()` cross-species ortholog mapping when `species` differs from `"Homo_sapiens"`.
  * Added `RunSCENICPlus()` for the SCENIC+ multi-omics workflow from Seurat objects, with a direct Python launcher, parallelized processing, and result readback.
  * Added `RunGRNBoost2()` and `RunGENIE3()` as standalone GRN modules wrapping the Arboreto Python implementations, with Seurat/matrix methods and `scenic_flt_adj()` target filtering shared with `RunSCENIC()`.
  * Added `RunCisTarget()` for standalone cisTarget motif enrichment and regulon construction from TF-target tables, with Python and R execution paths.
  * `PrepareDB()`: Added `data_dir` to parse locally downloaded single-file database sources (Broad MSigDB JSON, CSPA, Surfaceome, SPRomeDB, CORUM, JASPAR, ENCODE, TFLink, hTFtarget, TRRUST, CellTalk, CellChat) into the reusable `R.cache` database cache, avoiding repeated downloads.
  * `GeneSetScoring()`: Added C++ backends `zscore_dense()` and `plage_dense()` for Z-score and PLAGE gene-set scoring, plus PLAGE score orientation via z-score dot product for deterministic SVD signs. `ssgsea_rank_dense()` now accepts a `normalize` parameter. `aucell_auc_sparse()` gained a sparse `ctxcore` algorithm option.
  * `RunSCENIC()`: Added a `genome` parameter for automatic cisTarget reference selection while keeping `species` aligned with the package-wide Latin-name convention. Human references now support the default `"hg38"` v10 databases and `"hg19"` v9 databases, while the human TF list is cached with the genome-neutral name `allTFs_hgnc.txt`.
  * `RunSCENIC()`: Inlined the single-call `scenic_grn_inputs_changed` helper; shortened long internal function names (`scenic_dl_refs`, `scenic_flt_adj`, `scenic_def_mc`, `scenic_sel_mc_res`, `scenic_prep_gene_arg`) and synced cross-file calls in `RunGRN.R`.
  * `RunSCENICPlus()`: Inlined single-call helpers `scenicplus_read_optional_table`, `scenicplus_motif_tf_names`, `scenicplus_eregulon_table`; removed unused `scenicplus_eregulons`.
  * `FeatureHeatmap()`, `DynamicHeatmap()`, `GroupHeatmap()`: Compact long multi-term feature annotations from databases such as MSigDB and Reactome before drawing heatmap legends for cleaner display.
  * `CellCorHeatmap()`: Added `legend.position` parameter.
  * `GSVAPlot()`: Added `Database` column to enrichment results for consistent downstream filtering.
  * `RunMonocle2()`: Support custom root cells via `root_cells` parameter.
  * `RunDimsEstimate()`: Switched the default dimension-selection route to a scree-based ensemble of broken-stick, elbow, cumulative-variance, and marginal-gain criteria; the previous `intrinsicDimension` route remains available via `method = "intrinsic"` or can be combined with `method = "ensemble"`.
  * Added `srt_to_h5ad()` for writing Seurat objects to `.h5ad` files, complementing the existing AnnData-to-Seurat conversion utilities for scanpy interoperability.
  * Added `RunSciBet()` for native SciBet-style single-cell annotation from reference/query Seurat objects or expression matrices, with prediction metadata written back to Seurat.
  * Added `RunRareQ()` for RareQ rare-cell population detection from Seurat objects, including automatic neighbor construction through `DefaultReduction()`, metadata writeback, `CellDimPlot()` examples, and detailed result storage in `srt@tools[["RareQ"]]`.
  * Added `FerrisWheelPlot()` for up/down gene-count visualization from pathway enrichment results or pre-summarized count tables, with automatic `RunEnrichment()` result summarization, SCOP palette defaults, scalable outer donuts, title-cased pathway labels, and configurable text outlines.
  * Added `PalantirTrajectoryPlot()` for branch-aware Palantir trajectory visualization on Seurat embeddings, including pseudotime interval filtering, branch-probability path fitting, optional loess smoothing, branch-selection coloring, and layer-only return support.
  * Added `BranchStreamPlot()` for branch-aware pseudotime density ribbons from cell-state annotations and lineage pseudotime columns.
  * Added `RunMetaCell()` and `MetaCellPlot()` for metacell construction, original-cell to metacell mapping, metacell count output, and Seurat-compatible visualization.
  * Added `RunmcRigor()` for detecting dubious metacells or optimizing metacell partition granularity with runtime-optional `JSB-UCLA/mcRigor` installation, metadata writeback, and detailed result storage in `srt@tools[["mcRigor"]]`.
  * Added `RuntAge()` for tAge transcriptomic aging-clock prediction from Seurat pseudobulk, `ExpressionSet`, or matrix inputs, with runtime-optional `Gladyshev-Lab/tAge` preprocessing, R EN model caching from `mengxu98/datasets`, Python BR fallback support, and `tAgePlot()` visualization through `thisplot::StatPlot()`.
  * Added `RunscFEA()` for scFEA metabolic flux estimation from Seurat objects, with cached M168 resources from `mengxu98/datasets`, flux and balance assays, and `scFEAHeatmap()`, `scFEAVolcanoPlot()`, and `scFEABalanceBarPlot()` visualization helpers.
  * Added `RunAugur()` for Augur cell-type perturbation prioritization from Seurat objects, with an optimized backend and metadata/tool-slot writeback.
  * Added `RunCCC()` to run CellChat, CellphoneDB, and LIANA through one scheduler and rebuild a unified `srt@tools[["CCC"]]` bundle. These CCC wrappers now expose `backend = c("cpp", "r")` for scop post-processing and plotting-table aggregation while preserving each upstream package's inference logic.
  * Added `RunSCENIC()` for a SCENIC workflow from Seurat objects, including GRNBoost2/`scenic ctx` execution, regulon conversion, multi-core AUCell batch scoring, and storage of regulon activity scores as a Seurat assay plus detailed results in `@tools`. Metacell SCENIC workflows should now pass a metacell-level object from `RunMetaCell()` directly to `RunSCENIC()` instead of using an internal `group.by` aggregation path.
  * `RunSCENIC()` now uses a single `backend` argument to select the full SCENIC execution path. The public motif-fallback, AUCell backend, AUCell C++ strategy, and AUCell batch-size switches were removed from the SCENIC wrapper so GRNBoost2, cisTarget pruning, and AUCell scoring are controlled consistently by `backend = "cpp"` or `backend = "python"`.
  * Added `SCENICPlot()` to calculate regulon specificity scores from SCENIC activity and plot the top regulons for each metadata group.
  * `SCENICPlot()` heatmaps now expose `rss_scale`, `heatmap_limits`, and `heatmap_order`, allowing RSS and activity heatmaps to use matched row-wise z-score color scales and stable group-block row ordering when requested.
  * Added `RunScissor()` and `ScissorPlot()` for Scissor phenotype-associated cell selection from Seurat and bulk/SummarizedExperiment inputs, with an optimized backend, upstream-package backend, Seurat metadata/tool writeback, embedding plots, heatmaps, and status-composition summaries.
  * Added `RunscTenifoldNet()` for condition-level scTenifoldNet network comparison from matrices or Seurat groups using `cailab-tamu/scTenifoldNet`.
  * `RunscTenifoldKnk()` now uses the optimized path directly for QC, network-ensemble construction, tensor denoising, manifold alignment, and differential-regulation summaries; `scTenifoldKnkPlot()` includes QQ, effect-size, network, manifold, volcano, and upset result views. Added `scTenifoldNetPlot()` for condition-level scTenifoldNet QQ, effect-size, network, and manifold views.
  * Added `RunCIBERSORT()` and `ImmuneAbundancePlot()` for immune cell abundance deconvolution and visualization from bulk-like expression matrices, with C++ benchmarking support and lazy optional backend handling.
  * Added `NMFHeatmap()` for cell- or feature-level NMF similarity heatmaps with optional enrichment annotations and progress logging for large render jobs.
  * `RunPAGA()` now supports a C++ backend for the standard PAGA connectivity graph and uses it by default; `backend = "python"` remains available for RNA-velocity transitions, PAGA-initialized embeddings, plotting side effects, and DPT pseudotime.
  * `RunSCVELO()` now includes an optimized C++ backend for stochastic velocity embedding, compatible with `VelocityPlot()` and substantially faster than the equivalent R reference calculation. The Python backend remains the default for the full scVelo workflow.
  * `VelocityPlot()` raw, grid, and stream visualizations were validated with the `RunSCVELO(backend = "cpp")` velocity embeddings.
  * `PrepareEnv()` now supports explicit external wrapper modules (`"scmalignantfinder"`, `"secact"`, `"scpagwas"`, and `"external_wrappers"`) to preflight optional GitHub/Python dependencies without making those upstream tools hard package dependencies.
  * `PrepareEnv()` now supports `modules = "scenic"` as a standalone Python 3.10 environment (`scenic_env` by default) with SCENIC 0.12.1 and numpy 1.23.5, avoiding conflicts with the default `scop_env`.
  * Added `ConvertHomologs()` for homologous feature conversion in `Seurat`, `matrix`, and `Matrix` objects. The function uses `GeneConvert()` for arbitrary Ensembl/biomaRt-supported species pairs, collapses duplicated target homologs by summing expression values, preserves Seurat cell metadata and spatial images, and stores the mapping table in `@tools$ConvertHomologs`.
  * Added `RunCytoSPACE()`, an R/C++ implementation of the default CytoSPACE spot-level assignment workflow. The C++ backend uses spot-capacity graph construction and precomputed Pearson correlation matrices, stores detailed results in `srt@tools[["CytoSPACE"]]`, and writes summary metadata columns with the requested prefix.
  * Added `RunRCTD()` and `RunSPOTlight()` spatial deconvolution wrappers for estimating spot-level cell-type proportions from spatial Seurat objects and single-cell references, with standardized metadata columns compatible with `SpatialSpotPlot()` and `standard_scop(spatial_deconv_method = ...)`.
  * Added `RunSpatialEcoTyper()` as an optional SpatialEcoTyper wrapper for single-sample discovery, multi-sample conserved SE discovery, pretrained SE recovery, and SE abundance deconvolution. Results are written back to Seurat metadata and `srt@tools`, and `SpatialEcoTyperSpatialPlot()`/`SpatialEcoTyperCompositionPlot()` provide SCOP-styled visualization helpers.
  * Added `RunSpatialGradientFeatures()` and `SpatialGradientPlot()` for SPATA2-compatible spatial gradient feature screening, including a native C++ screening path, trajectory or annotation-based screening modes, model-fit summaries, and spatial/line/profile visualizations.
  * Added lightweight `semla` spatial wrappers: `RunSemlaSpatialNetwork()`, `RunSemlaLocalG()`, `RunSemlaRegionNeighbors()`, and `RunSemlaRadialDistance()`.
  * Added Giotto integration helpers for spatial workflows, including standalone result wrappers (`RunGiottoCluster()`, `RunGiottoSpatialGenes()`, `RunGiottoSpatialModules()`, `RunGiottoCellProximity()`), `GiottoPlot()` result visualizations, and the `SeuratToScopGiotto()` / `RunGiottoWorkflow()` object workflow for keeping Giotto results separate from Seurat until explicitly written back.
  * Added `SpatialSpotPlot()` for spatial visualization, including examples that show both tissue annotations and downstream CytoSPACE assignment results.
  * Added a shared C++ progress helper in `src/log_message.h` for long-running C++ loops. CytoSPACE assignment, scTenifold tensor decomposition, proportion permutation/bootstrap, and sample-level proportion bootstrap now report progress with the same timestamped format as `thisutils::log_message()`.

* **fix**:
  * `RunscPagwas()` now prefers the upstream `scPagwas_main2` runner, accepts validated GWAS data-frame or file inputs, forwards single-cell, cell-type, assay, and return controls explicitly, and applies local Seurat 5 compatibility without mutating optional-package namespaces. `RunscPaGWAS()` is retained as a deprecated alias, and `PlotScPagwas()` provides score and significance embedding plots.
  * `SpatialBackendStatus(api_check = TRUE)` now validates the actual namespace functions consumed by SpatialQM, SpotSweeper, CARD, SpatialEcoTyper, BayesSpace, smoothclust, MERINGUE, SPARK-X, nnSVG, and semla. CARD diagnostics resolve either CARD or CARDspa, including CARDspa's namespace-only constructor, and `RunCARD()` routes both legacy dotted and current underscore constructor arguments. Semla diagnostics select producer-specific symbols for its network, Local G, region-neighbor, and radial-distance wrappers.
  * Migrated `RunSpaNorm()`, `RunSpotSweeper()`, `RunRCTD()`, `RunCARD()`, `RunSpatialDWLS()`, `RunBANKSY()`, `RunCytoSPACE()`, `RunSmoothClust()`, `RunMERINGUE()`, `RunSpatialVariableFeatures()`, `RunSpatialGradientFeatures()`, `RunSpatialNeighborhood()`, `RunStatialKontextual()`, `RunSpatialIntegration()`, and `RunMistyR()` to raw-coordinate defaults. Their documentation now distinguishes coordinate-valued radii, bandwidths, and thresholds from unitless neighbor counts and preserves `"legacy_display"` as an explicit compatibility choice.
  * Spatial deconvolution registry entries now point RCTD, CARD, SPOTlight, and SpatialDWLS to the stored-result-only `SpatialDeconvolutionPlot()`, with strict spot identity validation and shared 0-1 proportion scales. C-SIDE is no longer misrepresented as a proportion result, and `RunRCTD()` now supports custom schema-v1 tool keys.
  * `RunSpatialIntegration(method = "SpatialMNN")` now targets the upstream package's actual `atlasClustering` namespace and runs `stage_1()` then `stage_2(method = "MNN")`, then strictly maps the returned sample-cluster table back to spots without invoking the backend's plotting side effects. PRECAST now consumes the SCOP-selected feature set, calls `SelectModel()` with its real argument name, and extracts its S4 domain/embedding results. BASS discovery now rejects the unrelated CRAN package with the same name and targets `zhengli09/BASS`.
  * `PrepareDB()`, `RunEnrichment()`, `RunGSEA()`, and `GeneConvert()` now accept Latin species names with spaces, dots, or hyphens (for example, `species = "Bos taurus"`) by normalizing them to the package's underscore form before database lookup and gene ID conversion.
  * `RunSCENIC()`: The C++ backend now writes the final filtered TF list to `<prefix>_regulators.txt` and records it in `@tools$SCENIC$files`, matching the persisted output bundle users expect beside the C++ adjacency, regulon list, and activity-score files.
  * `RunSCENIC()`: The C++ backend now labels positive regulons as `TF(+)` to match pySCENIC naming and supports optional negatively correlated regulons as `TF(-)` with `include_negative_regulons = TRUE`.
  * `SCENICPlot()`: `rss_heatmap` and `activity_heatmap` now return the drawable heatmap plot in `$plot` instead of the full `FeatureHeatmap()`/`GroupHeatmap()` result list, so direct PDF export draws the heatmap rather than printing metadata.
  * `SCENICPlot()`: `plot_type = "rss_rank"` now keeps all requested `top_n` and highlighted TF labels by default through `label_max_overlaps = Inf`; users can lower `label_max_overlaps` to let `ggrepel` drop crowded labels.
  * `RunCytoTRACE()`: Default to the official `CytoTRACE2::cytotrace2()` R package via `backend = "r"` and keep the native `scop` R/C++ implementation available through `backend = "cpp"`. The C++ backend remains independent of the official R package and reads model files from `PrepareDB(db = "CytoTRACE2")` by default. Added C++ subsample progress through `thisutils::parallelize_fun()`, in-memory model-data caching, dense BLAS multiplication for the near-dense CytoTRACE2 background matrices, C++ rank/log2-CPM preprocessing, top-k-only kNN smoothing, fewer exported Rcpp helper entrypoints, and official-R-vs-C++ parity/benchmark scripts covering mouse, human, explicit `data_dir`, multiple data sizes, and 1/2/4-core runs.
  * Cell-cell communication wrappers now preserve upstream result semantics more strictly: `RunNichenetr(mode = "aggregate_cluster_de")` passes receiver arguments with the current NicheNet API, `RunMultiNichenetr()` respects user-provided contrast tables while keeping original cell-type labels in stored parameters, and CellChat plots with condition-specific filters read from the original CellChat object instead of the cached unified CCC table.
  * `RunAugur()`: The `backend = "cpp"` path now performs Augur variance and random feature selection inside `scop`, so it no longer requires the GitHub-only `Augur` package unless `backend = "r"` is requested. C++ subsample results are accumulated before row binding to reduce repeated table growth, and failed C++ cell-type tasks now stop with a direct result-structure error instead of falling through to a secondary missing-`metric` error.
  * `RunPalantir()`: Fixed saved plot generation by passing `plot_format` through to Python, saving each embedding from its own matplotlib figure, and accepting scalar `early_group` values such as `"8"` when selecting the starting group.
  * `SCENICPlot()`: Explicit `features` in SCENIC heatmaps now keep the user-supplied regulon order, and `activity_heatmap` aligns `feature_split` to the resolved and displayed regulons.
  * `SCENICPlot()`: `plot_type = "activity_dim"` and `"activity_violin"` now respect all explicitly supplied `features`, instead of applying the six-regulon default preview limit.
  * `SCENICPlot()`: `plot_type = "target_bar"` now respects all explicitly supplied `features`, instead of applying the four-regulon default preview limit.
  * `RunscFEA()`: Configure Python thread and OpenMP runtime variables before dependency checks/imports to avoid torch runtime conflicts on macOS and mixed BLAS environments.
  * `scFEAHeatmap()`: Use the package heatmap palette conventions (`RdBu`, `Chinese`, and `simspec`), border/raster defaults, non-clustered group summaries, and `ht_params` overrides for closer consistency with `GroupHeatmap()`.
  * `GeneSetScoring()`: Route dense Gaussian and Poisson GSVA scoring through the C++ KDE path instead of immediately falling back to sparse-exact scoring, keeping the optimized z-score/rank workflow active for dense inputs.
  * `RunHarmony2()`: Added compatibility with Harmony 2.0 objects by directly trying both legacy fields (`Z_corr` / `R`) and callable methods (`getZcorr()` / `getR()`), including module methods that are not listed by `ls()`.
  * `srt_append()`: Align variable-feature metadata by feature name when appending into an existing Assay5 with a different feature universe, avoiding row-count replacement errors after integration workflows.
  * `EnrichmentPlot()` and `GSEAPlot()`: Resolve database aliases before applying `group_use`, and report selected groups with no enrichment rows directly instead of misreporting the database as missing.
  * `RunPalantir()`: Avoided a macOS/Apple Silicon crash caused by `umap-learn` importing TensorFlow's ParametricUMAP path during Palantir diffusion-map construction.
  * `RunCellTypist()`: Avoided a full AnnData-to-Seurat roundtrip for the common metadata-only annotation path; `RunCellphoneDB()` now avoids repeated Python environment checks and expands result tables with a vectorized path.
  * `RunCellphoneDB()`: Replaced the internal manual homolog-expression conversion path with `ConvertHomologs()`, keeping expression-object conversion behavior consistent across the package.
  * `GeneConvert()` examples now direct expression-object homolog conversion to `ConvertHomologs()` instead of showing manual `geneID_expand` aggregation.
  * `PrepareDB()`: Added compatibility with older `GOSemSim::godata()` signatures that use `OrgDb` instead of `annoDb`, avoiding GO semantic-data preparation failures in mixed Bioconductor environments.
  * `PrepareDB()`: Fixed direct MSigDB subcollection requests such as `db = "MSigDB_H"` or `db = "MSigDB_MH"` by resolving the base MSigDB species metadata before selecting the Broad release.
  * `PrepareDB()` and `AnnotateFeatures()`: Normalize legacy MSigDB caches whose feature column was stored as `symbol.ensembl_id`, and ensure ID-type conversion uses a single existing source ID column to avoid `switch()` errors when annotating MSigDB features by `symbol`.
  * `DEtestPlot()`, `VolcanoPlot()`, `DEtestManhattanPlot()`, and `DEtestRingPlot()`: Added `label.by` to choose automatic top-gene labels by adjusted p-value, p-value, detection-rate difference, or log2 fold change. Volcano and Manhattan plots now keep displayed positions and colors tied to the raw `avg_log2FC` values, and Manhattan plots default to no vertical jitter while allowing the centered group track size to be overridden with `group_track_width` and `group_track_height`.
  * `DynamicHeatmap()`, `FeatureHeatmap()`, and `GroupHeatmap()`: Compact long multi-term feature annotations from databases such as MSigDB and Reactome before drawing heatmap legends, and keep `DynamicHeatmap()` cluster annotations such as `RNA_snn_res.0.8` discrete even when stored as numeric metadata.
  * Cleaned up package-check issues by declaring missing namespace imports and aligning Rd argument documentation for recently updated wrappers.
  * Optional wrapper dependencies are checked at function entry with `check_r()` instead of silently skipping examples or adding unnecessary hard dependencies.

* **docs**:
  * Updated the pkgdown reference grouping for spatial analysis, spatial visualization, data conversion, and composition-analysis functions.
  * Updated the pkgdown reference index for CNV, ESTIMATE, SPOTlight, semla, Giotto, spatial-gradient, and immune-abundance workflows.
  * Updated `RunCytoSPACE()` examples to use real bundled data, convert mouse reference data with `ConvertHomologs()`, and visualize assignment results.
  * Refreshed examples to use bundled real package data and removed unnecessary `dontrun` wrappers from examples that do not require Python or external command-line tools.

* **data**:
  * Added `visium_human_pancreas_sub`, a Visium human pancreas spatial example dataset ([GSE254829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254829)) with spatial image data and CODA-derived annotations, for spatial analysis.
  * Added a Xenium package-data import script for preparing curated Xenium example resources outside the installed package payload.

# scop 0.8.9

* **fix**:
  * `RunGSVA()`: Removed redundant dense `as.matrix()` conversion in single-cell mode for `backend = "r"`; sparse expression matrices now stay sparse through row filtering, reducing peak memory and avoiding unnecessary dense materialization.
  * `RunCellQC()`: Replaced `Seurat::SplitObject()` with lazy cell-name-based subsetting, avoiding full-object duplication per split group. Added caching of `GetAssay()` gene names, `GetAssayData5()` counts, and per-cell QC metrics (`nCount` / `nFeature`) outside the species-check loop to eliminate repeated data extraction.
  * `RunDEtest()`: Pre-computed cell index mapping (`split(names(cell_group), cell_group)`) for paired-marker tests, replacing per-pair `which()` calls with O(1) list lookups.
  * `AnnotateFeatures()`: Replaced per-detail `sapply(… "[")` with type-stable `vapply(…, character(1))` to avoid implicit list-to-character coercion.
  * `run_scomm()`: Deferred dense conversion for reference and query subsetting by removing premature `as.matrix()` calls; sparse matrices are subset first, then densified only when required.
  * `RunUMAP2()`: Replaced `isSymmetric(as.matrix(graph))` with sparse-aware `Matrix::isSymmetric(graph)` in symmetry checks and graph subset sampling, avoiding dense materialization of large neighbor graphs.
  * `RunUMAP2()`: Replaced `apply(as.matrix(graph), 2, order)` in the uwot-predict path with the internal C++ sparse column top-k helper (`run_sparse_topk_by_column()`), avoiding full dense conversion and column-wise R-level `apply()`.
  * `RunDimsReduction()` (PCA centering): Uses `SeuratObject::LayerData(…, features = features)` to read only HVF rows from `scale.data` instead of loading the full matrix followed by manual subsetting.
  * `RunMDS()`: Removed `as.matrix()` before `Matrix::t()`; `proxyC::dist` now receives the sparse matrix directly, avoiding dense conversion of large count matrices.
  * `integration.R` (fastMNN): Removed three `as.matrix()` calls passed to `batchelor::fastMNN()`, which accepts sparse matrices directly via `SingleCellExperiment`.
  * `standard_scop()`: Replaced full `scale.data` matrix load (via `GetAssayData5()`) with direct layer access (`SeuratObject::GetAssayData` for Assay5; `@scale.data` for Assay) to retrieve only rownames when checking whether HVFs have been scaled, substantially reducing peak memory during the ScaleData decision step.
  * `GetAssayData5.Assay()`: Fixed parameter naming to use positional matching for the slot/layer argument, so `GetAssayData5()` correctly retrieves `counts`, `data`, and `scale.data` layers in both Seurat v4 (`slot`) and v5 (`layer`). Previously the named `layer` argument was silently dropped by SeuratObject v4, always returning the `data` slot regardless of the requested layer.
  * `CSS_integrate()`: Added `Assay5` guard before `SeuratObject::JoinLayers()` to avoid errors with Seurat v4 `Assay` objects, which do not support layered storage.
  * `RunDimsReduction()`: Added Seurat v4 fallback for PCA centering — `SeuratObject::LayerData(…, features = …)` is Assay5-only; v4 `Assay` objects now use `GetAssayData5()` with manual feature subsetting.
  * `integration_scop()`: Added early detection for v5-only integration methods (`CCA`, `RPCA`, `fastMNN5`, `Harmony5`, `scVI5`) to provide a clear, actionable error message when used with Seurat v4 `Assay` objects, rather than failing deep in the call stack.

* **feat**:
  * Added `loom_to_srt()` for pure-R loom-to-Seurat conversion via `rhdf5`, preserving velocity `spliced` and `unspliced` layers as assays without initializing Python, and added Python-backed `loom_to_adata()` for users who need AnnData output.
  * Added `RunBulk()` as a unified bulk-strategy entrypoint with method-vector selection for bulk DE, deconvolution, and cell-type-specific DE workflows.
  * Added method-specific bulk runners for `de_limma_voom`, `de_edgeR_qlf`, `de_DESeq2`, `de_dream`, `deconv_MuSiC`, `deconv_BisqueRNA`, `deconv_BayesPrism`, and `csde_TOAST`.
  * Standardized bulk results under `Bulk$results$de`, `Bulk$results$deconv`, and `Bulk$results$csde`, keeping DE outputs compatible with existing `DEtestPlot()`, `RunEnrichment()`, `RunGSEA()`, `GroupHeatmap()`, and `FeatureHeatmap()` data flows.
  * `RunBulk(run_enrichment = TRUE, run_gsea = TRUE)` now filters bulk DE rows with pathway thresholds and mirrors successful pathway results to `Enrichment_Bulk_wilcox` and `GSEA_Bulk_wilcox`, so `EnrichmentPlot()` and `GSEAPlot()` can read them through the standard `srt@tools` contract.
  * Deconvolution and CSDE bundles now record their computational `engine` in `details`; the current deconvolution runners use explicit `backend = "internal"` SCOP profile fitting and `csde_TOAST` uses explicit `backend = "limma_interaction"`, so package backends can be wired in later without pretending to call external packages.
  * Added explicit `Remotes` entries for `mengxu98/thisplot` and `mengxu98/thisutils` so source installs can resolve the minimum imported versions required by SCOP. Optional bulk engines such as `DESeq2` and `variancePartition` are installed on demand through `check_r()` and called lazily through `get_namespace_fun()`.
  * Added `method_args` to expose method-specific tuning parameters without expanding the public API surface.

# scop 0.8.8

* **deps**:
  * Updated the minimum dependency versions to `thisplot (>= 0.3.8)` and `thisutils (>= 0.4.5)`, and removed local copies of helpers now provided upstream (`thisplot::annotate_quadrants()`, `thisplot::clip_symmetric_range()`, `thisutils::collapse_sparse_rows()`).

* **docs**:
  * Updated pkgdown reference grouping for the cell-cycle workflow.

* **fix**:
  * Added optional wrappers for `RunDorothea()`, `RunBayesSpace()`, and experimental `RunscTenifoldKnk()`. `RunscTenifoldKnk()` keeps the upstream `scTenifoldNet` workflow but fixes the QC gene-filter assignment in the local path, uses an equivalent covariance/downdate path with direct sparse matrix construction, selection-based quantile thresholding, and controlled per-gene eigensolver parallelism for large `pcNet()` network construction, and uses C++ helpers for tensor decomposition, manifold matrix construction, directionality, and differential-regulation distance calculations. The C++ tensor-decomposition path now computes MTTKRP updates directly instead of materializing four dense unfolding matrices.
  * `PrepareEnv()`: Added explicit support for `mamba` and `micromamba` executables in addition to `conda`, including command-name resolution, automatic package-managed micromamba download when `conda = "micromamba"` is requested and the command is not on `PATH`, manager-aware logging, micromamba env path detection, and ToS handling that only runs for standard conda. Package-managed micromamba environments now use a no-space cache root to avoid `micromamba run -p` path parsing failures. `PrepareEnv()` now also stops early if reticulate has already initialized a different Python executable, avoiding unsafe in-session switches between conda and micromamba. The Python stack now also pins `setuptools < 81` so legacy packages such as `trimap` can still import `pkg_resources`.
  * `PrepareEnv()`: The default `modules = NULL` environment no longer includes `scomm`; install it explicitly with `modules = "scomm"` because the TensorFlow/scOMM stack conflicts with the default JAX/scVI stack through incompatible `ml-dtypes` requirements.
  * `CheckDataList()` / `standard_scop()`: Removed an unnecessary all-feature `ScaleData()` call immediately after `LogNormalize`. Downstream workflows still scale the selected HVFs before PCA, but large raw-count inputs no longer create a dense all-gene `scale.data` intermediate during checking.
  * `srt_reorder()`: Replaced per-cluster repeated sparse subsetting and `rowMeans()` with a single sparse group-membership matrix multiplication, speeding the cluster-reordering step used by `standard_scop()` while preserving average expression values.
  * `FindExpressedMarkers()`: Added a sparse fold-change and detection-rate prefilter path for common sparse RNA inputs, avoiding dense all-feature expression materialization before marker filtering while preserving default results. Supplied features are now intersected with the active layer before dense fallback scoring, which also keeps `scale.data` and partial-feature layers from passing absent rows into fold-change calculation.
  * C++-accelerated backends are now the default where available: `RunMetabolism()`, `RunGSVA()`, `CellScoring()`, `RunDynamicEnrichment()`, and `RunEnrichment()` now prefer `backend = "cpp"` while retaining `backend = "r"` for exact legacy/package behavior. `RunPermutation()` now uses its validated C++ implementation directly.
  * `RunMetabolism()` and `RunGSVA()`: Added a reusable C++ gene-set scoring backend for `ssGSEA`, `zscore`, and `plage`; `RunGSVA()` can now use `backend = "cpp"` for `method = "ssgsea"`, `method = "zscore"`, `method = "plage"`, and Gaussian- or Poisson-kernel `method = "gsva"`. PLAGE scores are oriented by the gene set mean z-score to avoid backend-dependent sign flips.
  * `RunMetabolism()` and `RunGSVA()`: Added `cpp_chunk_size` for the C++ GSVA kernel paths to reduce peak dense intermediate memory on large cell counts; `NULL` now auto-selects a chunk size for large matrices.
  * `RunEnrichment()`: Added an experimental `backend = "cpp"` ORA path using a C++ hypergeometric implementation for faster enrichment tables while keeping `backend = "r"` available as the clusterProfiler-compatible path.
  * `CellScoring()`: Added experimental `backend = "cpp"` support for Seurat module scoring by keeping control-gene sampling in R and moving sparse mean calculations to C++ code.
  * `RunPermutation()`: Uses the C++ permutation and bootstrap loops directly after validating they match the legacy R calculation for observed fractions and log2 fold-differences while running substantially faster.
  * `RunUMAP2()`: Added an internal C++ sparse column top-k helper for Graph inputs to speed extraction of precomputed neighbor indices/connectivities before calling `uwot`.
  * `RunKNNMap()` and `RunKNNPredict()`: Added an internal C++ dense column top-k helper for the raw KNN fallback to speed nearest-neighbor extraction from precomputed distance matrices.
  * `PseudotimeProjectionPlot()`: Reused the internal dense top-k helper for pseudotime KNN/gradient neighbor extraction from distance matrices.
  * Mapping/integration metrics: Added an internal C++ contingency-table helper for `metric_accuracy()`, `metric_macro_f1()`, `metric_purity()`, `metric_nmi()`, `metric_ari()`, `metric_weighted_recall()`, and `collect_mapping_metrics()`.
  * `RunMetabolism()`: Added a C++ backend for `method = "GSVA"` that keeps the Poisson-kernel GSVA scoring path but accelerates repeated count-kernel calculations.
  * `RunMetabolism()` and `CellScoring()`: Added an experimental C++ backend for AUCell scoring selected only through `backend = "cpp"`, while keeping the original R/AUCell implementation available via `backend = "r"` for exact package output. `RunDynamicEnrichment()` now passes this backend through to AUCell-based scoring.
  * `CellScoring(method = "AUCell")`: Fixed score-name assignment when only one feature list is scored, avoiding vector dropping before metadata column names are applied.
  * `RunMetabolism(method = "VISION")`: Fixed signature construction so in-memory gene sets are passed as `VISION::Signature` objects instead of being interpreted as file paths.
  * `DynamicHeatmap()` / custom enrichment workflows: Fixed incorrect term label display when user-provided `TERM2GENE` / `TERM2NAME` use non-standard column names such as `term` / `name`. Custom term annotations are now normalized consistently so that heatmap term labels show readable term names instead of fallback IDs (e.g. GO IDs). Related issue #160 (@Pineapple-wen6, @mengxu98).
  * Refactored repeated custom database input handling by introducing shared internal helpers for normalizing and assembling user-provided `TERM2GENE` / `TERM2NAME`, and applied them across `RunEnrichment()`, `RunDynamicEnrichment()`, `RunGSEA()`, `RunGSVA()`, and `PrepareDB()` to keep behavior consistent.
  * `integration_scop()`: Improved `ChromatinAssay` handling by standardizing ATAC reductions after integration, avoiding non-interactive small-batch blocking, clipping `TFIDF/rlsi` anchor dims to available cells, auto-switching `Harmony5` to legacy `Harmony`, and rejecting unsupported `Seurat` / `RPCA` ATAC paths with explicit messages.

# scop 0.8.7

* **feat**:
  * Added `RunDimsEstimate()` and `DimsEstimatePlot()` for intrinsic dimensionality estimation from reductions, integrated into `RunDimsReduction()` and `standard_scop()` to automatically select useful dimensions when `linear_reduction_dims_use = NULL`.
  * Renamed `RunDimReduction()` to `RunDimsReduction()` and updated downstream callers/documentation accordingly.
  * `standard_scop()`: When `linear_reduction_dims_use = NULL`, now uses estimated dimensions stored in the reduction (via `RunDimsEstimate()`) when available, with fallback to the first 50 dimensions.
  * Added Seurat v5 integration methods: `CCA_integrate()`, `RPCA_integrate()`, `fastMNN5_integrate()`, `Harmony5_integrate()`, and `scVI5_integrate()` via `Seurat::IntegrateLayers()`, and exposed them through `integration_scop()`.
  * Added `Coralysis_integrate()` and exposed Coralysis through `integration_scop()`, following the official Seurat v5 compatible workflow via `SingleCellExperiment`.
  * Added `h5ad_to_srt()` for reading `.h5ad` files directly into `Seurat` objects via `scanpy.read_h5ad()`, with automatic CSR/float64 coercion to avoid reticulate conversion issues. Layers that fail conversion are gracefully skipped and reported.
  * `adata_to_srt()`: Improved robustness of layer conversion — each layer is now wrapped in `tryCatch()` so that individual failures no longer abort the entire conversion; skipped layers are reported as warnings.
  * `RunCellChat()`: Enhanced to support condition-specific analyses and pairwise merged comparisons via `group_column` and `group_cmp` parameters.
  * Added `RunCellphoneDB()` for running CellphoneDB cell-cell communication analysis on a `Seurat` object through the official Python package, with support for species conversion and results stored in `srt@tools[["CellphoneDB"]]`.
  * Added `RunNichenetr()` and `RunMultiNichenetr()` for running NicheNet and MultiNicheNet analysis on `Seurat` objects with standardized result storage.
  * Added unified cell-cell communication plotting functions `CCCStatPlot()`, `CCCHeatmap()`, and `CCCNetworkPlot()` for CellChat, CellphoneDB, NicheNet, and MultiNicheNet results.
  * Remove `CellChatPlot()`.

* **fix**:
  * `RunUMAP2()`: Fixed reduction lookup to check existing reduction names before falling back to `DefaultReduction()`, avoiding errors when the exact reduction name is already present.

# scop 0.8.6

* **feat**:
  * Added `RunGSVA()` for gene set variation analysis and `GSVAPlot()` for visualization of [GSVA](https://github.com/rcastelo/GSVA) results. Related issue #146 (@mengxu98).
  * `RunMetabolism()`
    * Added `RunMetabolism()` for single-cell metabolism pathway scoring with support for `AUCell`, `GSVA`, `ssGSEA`, and optional `VISION`, and `MetabolismPlot()` for visualization
    * `RunMetabolism()` now uses [scMetabolism](https://github.com/wu-yc/scMetabolism) KEGG / Reactome metabolism pathway definitions to identify metabolism-related terms, then rebuilds pathway gene sets from updated `PrepareDB()` annotations, enabling cached database updates and cross-species metabolism scoring.
    Related issue #146 (@mengxu98).
  * Added `RunDecontX()` and integrated optional `decontX` ambient RNA decontamination into `RunCellQC()`, including contamination metadata, optional decontaminated assay output, and threshold-based QC filtering. Related issue #147 (@mengxu98).

* **fix**:
  * `FeatureStatPlot()`: Fixed duplicated X/Y axis titles and main title when `stack = TRUE` with `theme_args` or `title`. Stack assembly now strips axis and plot titles from non-edge panels and draws a single shared ylab/title via gtable; theme styling from `theme_args` (e.g. `axis.title.y`, `title`/`plot.title`) is applied to these shared grobs. Related issue #145 (@PanSX-Dr).

* **docs**:
  * Updated the README badge/logo.

* **data**:
  * Removed the example datasets `ifnb_sub`, `ref_scHCL`, and `ref_scZCL`.

# scop 0.8.5

* **feat**:
  * Default palette changed from `"Paired"` to `"Chinese"`. See `thisplot::ChineseColors()` for details.
  * Unified the `cores` parameter across multiple functions (`DynamicHeatmap()`, `FeatureHeatmap()`, `GroupHeatmap()`, and `heatmap_enrichment()`), ensuring `cores` is correctly threaded through to `RunEnrichment()`.
  * `RunMonocle2()` and `RunMonocle3()`: Added `xlab` and `ylab` parameters, passed through to internal `CellDimPlot()` and `FeatureDimPlot()` calls.
  * `ListDB()`: Now supports multiple species in the `species` parameter simultaneously and adds `Species` and `DB` columns to the output data frame for clearer identification.
  * Moved `is_outlier` to `thisutils::is_outlier()`.
  * `configure_apple_silicon_env()` (*Python* function): Added OpenMP compatibility handling on macOS arm64 by prepending environment `lib` paths to `DYLD_FALLBACK_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` and preloading `libomp.dylib`/`libiomp5.dylib` via `ctypes.CDLL(..., RTLD_GLOBAL)` to reduce `scanpy`/`python-igraph` import conflicts.
  * `scVI_integrate()`: Unified parameter naming by renaming `num_threads` to `cores` for consistency across integration functions.
  * Added `find_neighbors_and_clusters()` and `run_nonlinear_reduction()` helper functions to unify operations across integration functions and reduce redundant code.

* **fix**:
  * `PrepareDB()`: 
    * TF database now uses [AnimalTFDB4](https://github.com/mengxu98/datasets/tree/main/AnimalTFDB4) as the data source.
    * MP database - the Web Archive URL year is now dynamically determined from the current date, and the archive snapshot date for all MP-related file downloads (`VOC_MammalianPhenotype.rpt`, `MGI_Gene_Model_Coord.rpt`, `MGI_GenePheno.rpt`) is extracted from the server index page to ensure consistent versioning.
    * hTFtarget data download URL has been changed from `"http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"` to `"https://guolab.wchscu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"`.
    * CSPA data download URL has been changed from `"https://wlab.ethz.ch/cspa/data/S1_File.xlsx"` to `"https://raw.githubusercontent.com/mengxu98/CSPA/main/S1_File.xlsx"`.

    Related issus #76 (@hwa2Hu), #139 (@pengding774-dot), #140 (@mengxu98).
  * `LIGER_integrate()`
    * Migrated to the `rliger` 2.x workflow (`rliger::runIntegration()` + `rliger::quantileNorm()` on `Seurat` object) and now prepares/uses the `ligerScaleData` layer via `rliger::scaleNotCenter()` before integration.
    * Updated argument naming from `LIGER_dims_use` to `liger_dims_use`, and removed legacy quantile-normalization parameter compatibility mapping (`ref_dataset`), keeping `reference` as the supported interface.
  * Optimized the installation of some *Python* packages on Apple Silicon devices.

* **docs**:
  * Updated README example code and visualizations. After regenerating figures, the package size was reduced by ~3 MB.

# scop 0.8.4

* **fix**:
  * `FeatureStatPlot()` / `ExpressionStatPlot()`: Fixed box and violin x-axis misalignment when `add_box = TRUE` with `split.by`. Groups with fewer than 2 observations are now filtered before violin density estimation (with a warning), and the violin layer uses a consistent `position_dodge(width = 0.9)` to match the boxplot. Related issue #123 (@oranges7). 

# scop 0.8.3

* **feat**:
  * `RunDynamicFeatures()`: Added [PreTSA](https://github.com/haotian-zhuang/PreTSA/) method for dynamic feature fitting. The [PreTSA](https://github.com/haotian-zhuang/PreTSA/) algorithm in `scop` has been re implemented to support parallelization for higher performance. Original research: [PreTSA: computationally efficient modeling of temporal and spatial gene expression patterns](https://doi.org/10.1186/s13059-026-03994-3). Use `fit_method = "pretsa"` for B-spline-based piecewise truncated spline analysis; `fit_method = "gam"` (default) keeps generalized additive models. PreTSA supports `knot` (`0` or `"auto"`) and `max_knot_allowed` when `knot = "auto"`. Relate issue #133.
  * `CellDimPlot()` and `FeatureDimPlot()`: Added `legend.title` parameter (default `NULL`) to control the legend title. When `NULL`, default titles are used (e.g. group name for `CellDimPlot`, feature or empty for `FeatureDimPlot`).
  * `ExpressionStatPlot()` and `FeatureStatPlot()`: Added `legend.title` parameter (default `NULL`) for single-legend plots. When `NULL`, the default title (e.g. `keynm` or feature/group name) is used. `FeatureStatPlot()` forwards `legend.title` to `ExpressionStatPlot()`.
  * Moved `StatPlot` function to `thisplot::StatPlot`.

* **fix**:
  * Differential expression plots (`VolcanoPlot()`, `DEtestManhattanPlot()`, `DEtestRingPlot()`): added `only.pos = TRUE` for positive-only visualization. `DEtestManhattanPlot()` now keeps the cell-type track centered at y = 0 and sizes the track from the nearest point distance around zero.
  * `DynamicHeatmap()` / `heatmap_enrichment()`: Fixed incorrect `db` handling when using custom `TERM2GENE`/`TERM2NAME`. Enrichment results with `Database = "custom"` could be incorrectly filtered by default `db` values (e.g. `"GO_BP"`), causing false "No term enriched using the threshold" warnings even when enrichment succeeded. Relate issue #133 (@1228849000).

# scop 0.8.2

* **feat**:
  * Add the `DEtestPlot()` function, which calls the original `VolcanoPlot()` and adds two plot types, Manhattan and Ring, controlled by the `plot_type` parameter (`c("volcano", "manhattan", "ring")`, default `"volcano"`). Add standalone functions `DEtestManhattanPlot()` and `DEtestRingPlot()` for direct use. Relate issue #121 (@ericavalentini).
  * Differential expression visualization (`DEtestPlot()`, `VolcanoPlot()`, `DEtestManhattanPlot()`, `DEtestRingPlot()`): added `res` parameter to accept existing DE results (data.frame). When `res` is provided, `srt` is ignored. Data processing supports: (1) `group1` or `cluster` column for grouped plots; (2) no grouping column for a single panel; (3) gene names from row names when `gene` column is missing. Relate issue #129 (@mengxu98).
  * `PrepareEnv()`: Update `version` parameter to specify the *Python* version of the conda environment. Default is `"3.11-1"` on *Windows* and `"3.10-1"` on *macOS* and *Unix*. Relate issue #103 (@PanSX-Dr).
  * `RunNMF()`: Add the `cores` parameter for `RunNMF()` and optimize the printed message.
  * Removed the setting in *Python* functions that prevents drawing functions from causing *R* crashes.

* **fix**:
  * `RunPalantir()`: Fixed unused `plot_format` parameter error. The parameter is now properly excluded from arguments passed to Python functions. Relate issue #126 (@christinejay990202-dev).

# scop 0.8.1

* **fix**:
  * `FeatureHeatmap()`: Fixed `group_palcolor` when a named vector is passed: the function previously used only the first color for all groups because `group_palcolor[[1]]` on a vector returns a single element. Now when `group.by` has length 1, a vector is automatically wrapped as a list; when `within_groups = TRUE`, `group_palcolor` is expanded in line with `group_palette`.
  * `GroupHeatmap()`: Same `group_palcolor` fix as `FeatureHeatmap()`: support for named-vector input and correct expansion when `within_groups = TRUE`.

* **docs**:
  * `group_by` modified to `group.by` to unify documentation, involving functions: `RunCellRank()`, `RunDEtest()`, `VolcanoPlot()`, `RunEnrichment()`, `EnrichmentPlot()`, `RunGSEA()`, `GSEAPlot()`, `RunPalantir()`, `RunPAGA()`, `RunWOT()`, `RunSCVELO()`, releated to [#120](https://github.com/mengxu98/scop/issues/120).

# scop 0.8.0

* **feat**:
  * `RunMonocle2()`: New function for performing [Monocle2](https://github.com/mengxu98/monocle) trajectory analysis with support for various dimensionality reduction methods (DDRTree, ICA, tSNE, SimplePPT, L1-graph, SGL-tree). Uses the fixed version of monocle2 from [mengxu98/monocle](https://github.com/mengxu98/monocle).
  * `RunMonocle3()`: New function for performing [Monocle3](https://github.com/cole-trapnell-lab/monocle3) trajectory analysis with support for cell ordering, trajectory learning, and pseudotime computation.
  * `RunCytoTRACE()`: New `scop` implementation for running CytoTRACE 2 analysis to predict cellular potency scores and categories (Differentiated, Unipotent, Oligopotent, Multipotent, Pluripotent, Totipotent) with support for human and mouse species.
  * `CytoTRACEPlot()`: New function for visualizing CytoTRACE 2 analysis results.

# scop 0.7.9

* **docs**:
  * Optimized `@inheritParams` usage to reduce redundant parameter definitions.
  * Unified parameter naming: renamed `group_by` to `group.by` for consistency with Seurat conventions in `RunMonocle3()`, `RunPalantir()`, `RunSCVELO()`, `RunWOT()`, `CellCorHeatmap()`, and `RunDynamicFeatures()`.

# scop 0.7.8

* **fix**:
  * `RunPalantir()`: Fixed unused `plot_format` parameter error. The parameter is now properly excluded from arguments passed to Python functions. Relate issue #114 (@Moonerss).
  * `RunSCVELO()`: Fixed PAGA computation error by replacing `scv.tl.paga` with `sc.tl.paga` (scanpy implementation) for better stability. The function now uses the same PAGA implementation as `RunPAGA()` function.
  * `RunCellRank()`: Fixed GPCCA Schur decomposition error by adding fallback mechanism. When `brandts` method fails with "subspace_angles" error, the function automatically tries `krylov` method. If both methods fail, it automatically switches to CFLARE estimator for more robust computation.
  * `RunCellRank()`: Fixed `recover_dynamics` error by ensuring `velocity_graph` and `velocity_graph_neg` are properly set before calling `scv.tl.recover_dynamics()` for latent time computation.

* **feat**:
  * `adata_to_srt()`: Removed automatic removal of "X_" prefix from dimensionality reduction names in `obsm` keys. The function now preserves original reduction names as they are stored in AnnData objects.

* **data**:
  * Reducing the size of `pancreas_sub` example dataset.

# scop 0.7.7

* **feat**:
  * `adata_to_srt()`: Enhanced to support multiple AnnData object types including *Python* AnnData objects (from scanpy/reticulate), R6 AnnData objects from the `anndata` package (AnnDataR6), and R6 AnnData objects from the `anndataR` package (InMemoryAnnData). Added internal helper functions `get_adata_element()` and `get_adata_names()` for better compatibility. Relate issues #67 (@lisch7), #91 (@mengxu98) and [commit91#issuecomment](https://github.com/mengxu98/scop/issues/91#issuecomment-3659404993).

* **fix**:
  * `RunDEtest()`: Fixed error when comparing one cluster against multiple clusters using `group1` and `group2` parameters. Relate issue #111 (@zhaoxiaoyan9225).
  * `AnnotateFeatures()`: Fixed bug where the function would fail when processing GTF file annotations due to column name matching issues during data naming. The function now correctly handles column name intersections when merging annotation data.

# scop 0.7.6

* **feat**:
  * `RunDM()`: Added automatic PCA-based dimensionality reduction when using many features (>1000) to speed up diffusion map computation. The `npcs` parameter can be used to control the number of principal components used for pre-processing.

# scop 0.7.5

* **fix**:
  * `CellScoring()`: Fixed bug where the function failed to build results. Relate issue #98 (@SuperrNaruto).

* **feat**:
  * `RunDEtest()`: Fixed compatibility issue with [SeuratObject](https://satijalab.github.io/seurat-object/) 5.0.0+ by replacing deprecated `Assays()` `slot` argument with `LayerData()`. Relate issue #100 (@mattizecos).
  * `RunDM()`: Added automatic PCA-based dimensionality reduction when using many features (>1000) to speed up diffusion map computation. The `npcs` parameter can be used to control the number of principal components used for pre-processing.

# scop 0.7.3

* **feat**:
  * `RunCellTypist()`: New function for cell type annotation using the CellTypist method.
  * `CellTypistModels()`: New function for downloading and managing CellTypist pre-trained models.

# scop 0.7.2

* **feat**:
  * `RunCellRank()`: Performance optimizations and code improvements.

# scop 0.7.1

* **fix**:
  * `CellDimPlot()`: Fixed issue where NA values appeared in labels. Relate issue #93 (@12345nkjil).

# scop 0.7.0

* **feat**:
  * `PrepareEnv()`: Integrated `uv` as the primary *Python* package installer for improved installation speed.
  * `check_python()`: Now uses `uv` as the primary installation tool with `pip` as fallback, significantly improving package installation speed.
  * Added `find_uv()` and `install_uv()` internal functions for managing `uv` package manager installation and detection.

# scop 0.6.6

* **docs**:
  * Unified documentation format across all R functions:
    * Standardized return value tags: Changed all `@returns` to `@return` for consistency.
    * Unified parameter documentation: Replaced all `\code{value}` with Markdown backticks `` `value` `` format.
    * Standardized default value descriptions.
    * Added `@md` tags: Added `@md` tags to all functions using Markdown syntax in documentation.
    * Enhanced cross-references: Added `@seealso` links to related functions where appropriate.

# scop 0.6.5

* **feat**:
  * `PrepareEnv()`:
    * Added comprehensive environment variable configuration to prevent crashes when calling *Python* functions, including setting thread limits for OMP, OPENBLAS, MKL, NUMBA, and other libraries. This improves stability on all platforms, especially Apple silicon Macs.
    * Added `accept_conda_tos()` function to automatically accept conda Terms of Service for required channels, improving the conda environment setup process.
    * Fixed conda Terms of Service acceptance issue in `PrepareEnv()`. The function now automatically accepts conda Terms of Service for required channels, eliminating the need for manual acceptance. This addresses the issue reported in [#85](https://github.com/mengxu98/scop/issues/85).
  * Multiple Python-based functions (`RunPAGA`, `RunSCVELO`, `RunPalantir`, `RunCellRank`, `RunWOT`, `RunPHATE`, `RunPaCMAP`, `RunTriMap`): Enhanced message formatting and code improvements.
  * `PrepareSCExplorer()`: Fixed package version dependency issues with `shiny` and `bslib` compatibility. The function now properly handles `bslib` theme configuration to work with both `shiny` 1.6.0 and 1.7.0+, addressing compatibility errors reported in [#87](https://github.com/mengxu98/scop/issues/87).

* **fix**:
  * Improved code formatting and consistency across multiple functions.
  * Enhanced Python functions in `inst/python/functions.py` with better error handling and message formatting.

* **docs**:
  * Updated documentation for multiple functions to reflect code improvements.

# scop 0.6.2

* **feat**:
  * `CellChatPlot()`: Adjusted the size of saved figures for better file size optimization.

* **docs**:
  * Updated README.md to remove references to Monocle2 and Monocle3 (deprecated functions).

# scop 0.6.1

* **feat**:
  * `PrepareEnv()`: Improved message formatting and simplified log output for better user experience.
  * Added `get_conda_envs_dir()` helper function to centralize conda environment directory retrieval.
  * `integration_scop()`: Enhanced `integration_method` parameter definition with explicit method list for better code clarity.

* **fix**:
  * Moved `exist_python_pkgs()` function to `check_package.R` for better code organization.
  * Replaced direct `conda_info()$envs_dirs[1]` calls with `get_conda_envs_dir()` helper function for consistency.
  * `RunSCExplorer()`: Updated to use `thisplot::palette_list` and `thisplot::slim_data()` instead of `scop::palette_list` and `scop::slim_data()`.
  * Added `thisplot` to dependency checks in `RunSCExplorer()`.

* **docs**:
  * Updated documentation across multiple functions.

# scop 0.6.0

* **feat**:
  * `PrepareEnv()`: Enhanced with environment caching mechanism to avoid redundant environment preparation. Improved message formatting and error handling.
  * Python-based functions (`RunPAGA()`, `RunSCVELO()`, `RunPalantir()`, `RunCellRank()`, `RunWOT()`) now automatically call `PrepareEnv()` internally, eliminating the need for users to manually prepare the Python environment before using these functions.
  * `cluster_within_group2()`: New function for clustering within groups.
  * Multiple plotting functions: Replaced `geom_sankey()` with `ggsankey::geom_sankey()` for better Sankey diagram support.
  * Multiple functions: Replaced `:::` operator with `get_namespace_fun()` for safer namespace access.

* **fix**:
  * Removed `RunMonocle()` function and related documentation (`RunMonocle2.Rd`, `RunMonocle3.Rd`).
  * Removed `projection_functions.R` file (functions moved to other locations).
  * Replaced custom theme functions with `thisplot::theme_this()` (exported as `theme_scop()`).
  * Replaced direct `log_message()` calls with `thisutils::log_message()` for consistency.
  * Removed `palette_list` data object.

* **deps**:
  * Moved `cli` from `Suggests` to `Imports` for better message formatting support.
  * Added `thisplot` to `Imports` for theme and utility functions.
  * Added `ggsankey` to `Suggests` for Sankey diagram support.
  * Added remote dependencies: `theislab/destiny`, `mengxu98/thisplot`.
  * Removed unused dependencies: `Biobase`, `BiocGenerics`, `concaveman`, `DDRTree`, `glmGamPoi`, `hexbin`, `monocle`, `png`, `ragg`, `tidyr`.

* **docs**:
  * Updated documentation across multiple functions to reflect code refactoring.
  * Improved code organization and maintainability.

# scop 0.5.5

* **fix**:
  * Fixed `VelocityPlot()` function error in `plot_type = "grid"` mode: replaced vectorized arrow length with fixed-length arrows (using mean length) to resolve `vapply()` error that occurred when `grid::arrow()` received a vector instead of a single value, see [#72](https://github.com/mengxu98/scop/issues/72), [#74](https://github.com/mengxu98/scop/issues/74).

# scop 0.5.4

* **fix**:
  * Fixed parameter name error in `CheckDataType()` function calls: changed `data` parameter to `object` in `RunKNNMap()`, `RunSymphonyMap()`, `RunScmap()`, `RunPCAMap()`, and `RunSingleR()` functions, see [#68](https://github.com/mengxu98/scop/issues/68).
  * Fixed `SingleCellExperiment` object creation in `RunScmap()` and `RunSingleR()` functions: changed from coercing `SummarizedExperiment` to directly constructing `SingleCellExperiment` objects.

# scop 0.5.3

* **feat**:
  * `PrepareDB()`: Changed default `Ensembl_version` parameter from `103` to `NULL` for more flexible version handling.
  * Added *Python* version `log_message()` for Python-based functions (`RunSCVELO()`, `RunPAGA()`, `RunPalantir()`, `RunCellRank()`, `RunWOT()`) and added `verbose` parameter inheritance and improved message formatting through cli formatting.

* **fix**:
  * Delete `harmonizomeapi.py` file.
  * Move `scop_analysis.py` into a single `functions.py` file in `inst/python/` for better code organization and maintainability.

* **docs**:
  * Improved parameter documentation consistency.

# scop 0.5.1

* **docs**:
  * Improved reference formatting and consistency across multiple functions.
  * Enhanced documentation clarity and readability.

# scop 0.5.0

* **feat**:
  * `RunCellChat()`: New function to perform CellChat analysis for investigating cell-to-cell communication with support for human, mouse, and zebrafish species.
  * `CellChatPlot()`: New function to visualize CellChat analysis results with various plot types and customization options.
  * Multiple integration functions: Improved error messages and message formatting for better user experience.

* **deps**:
  * Added [CellChat](https://github.com/jinworks/CellChat) package dependency with remote repository `jinworks/CellChat`.

* **docs**:
  * Updated README.md with improved code formatting and examples.
  * Enhanced documentation for cell communication analysis functions.
  * Improved error messages and user guidance across integration functions.

* **fix**:
  * Removed some example figures to optimize package installation size.

# scop 0.4.0

* **feat**:
  * `RunProportionTest()`: New function to perform Monte-carlo permutation test for quantifying cell proportion differences between conditions.
  * `ProportionTestPlot()`: New function to generate proportion test plots with customizable significance thresholds and visualization options.
  * Multiple *Python*-based functions: add `\dontrun{}` blocks for Github workfolw checking.

* **docs**:
  * Added comprehensive documentation for new proportion testing functions.
  * Enhanced example documentation across multiple functions.
  * Updated package documentation and examples.

# scop 0.3.4

* **docs**:
  * Updated workflow examples and function documentation.

# scop 0.3.3

* **feat**:
  * Multiple functions: Improved parameter documentation formatting and consistency across the package.

# scop 0.3.2

* **feat**:
  * `GetFeaturesData()` and `AddFeaturesData()`: Enhanced argument clarity, added input validation, and standardized return values for `Seurat`, `Assay`, and `Assay5` objects.
  * `CellCorHeatmap()`: 
    - Renamed parameters: `query_cell_annotation` → `query_annotation`, `ref_cell_annotation` → `ref_annotation`.
    - Improved error message formatting through cli formatting.
    - Simplified variable assignments and improved readability.

* **docs**:
  * Comprehensive documentation updates across multiple functions including `AnnotateFeatures`, `CellDimPlot`, `CellStatPlot`, `FeatureStatPlot`, `GroupHeatmap`, `RunCellQC`, and others.
  * Improved parameter descriptions and function clarity.

# scop 0.3.1

* **feat**:
  * `EnrichmentPlot()` and `GSEAPlot()`: Removed conditional font face styling (`face = ifelse()` logic) for better text rendering consistency. Set the default value of `lineheight` from `0.5` to `0.7`.
  * Updated `check_r()` function for improved package checking functionality.
  * Updated reexports functionality.

* **docs**:
  * Updated documentation formatting and consistency.

# scop 0.3.0

* **feat**:
  * Fixed `segmentation faults` and `R crashes` on *M-series* MacBook when running *Python* functions.
  * `RunPAGA()`: Enhanced with *M-series* MacBook detection and automatic environment configuration.
  * `RunSCVELO()`: Added ARM64-specific optimizations to prevent crashes and ensure stable execution.
  * `RunCellRank()`: Implemented *M-series* compatibility with proper NUMBA configuration.
  * `RunPalantir()`: Added ARM64 support with single-threaded execution mode.
  * `RunWOT()`: Enhanced with *M-series* MacBook environment variable settings.
  * `RunTriMap()`: Added *M-series* MacBook compatibility for dimensionality reduction.
  * `RunPaCMAP()`: Implemented ARM64-specific environment configuration.
  * `RunPHATE()`: Added *M-series* MacBook support for non-linear dimensionality reduction.
  * `RunCellQC()`: Enhanced both scrublet and doubletdetection functions with ARM64 compatibility.

* **docs**:
  * Updated function documentation to reflect *M-series* MacBook compatibility.
  * Added technical notes about ARM64 architecture considerations.

# scop 0.2.9

* **fix**:
  * Fix bug for `RunSCVELO()`.

* **docs**:
  * Updated documentation for some functions.

# scop 0.2.7

* **feat**:
  * Added an internal function `.check_pkg_status()` to check if an *R* package is installed.
  * Update function `CheckDataType()` to *S4* class function.
  * Update function `standard_scop()`, make it more efficient.

* **data**:
  * Delete `lifemap` data, including: `lifemap_cell`, `lifemap_compartment` and `lifemap_organ`.
  * Reconstructed the sample data for both `panc8_sub` and `pancreas_sub`, retaining only the basic `Seurat` object.

# scop 0.2.6

* **feat**:
  * Added `remove_r()` function for easy remove *R* packages.
  * Rename function: `RemovePackages()` to `remove_python()`.
  * Removed other methods of installing *R* packages from the `check_r()` function, only retaining [pak::pak](https://pak.r-lib.org/reference/pak.html). 
  * Delete useless import packages: `BBmisc`, `BiocManager`, `covr`, `devtools`, `promises` and `withr`.
  * Optimize the structure of `_pkgdown.yml` file.

* **docs**:
  * Updated documentation for some functions.

# scop 0.2.5

* **feat**:
  * Rename function: `palette_scop()` to `palette_colors()`.

# scop 0.2.4

* **feat**:
  * Rename functions: `check_srt_merge()` to `CheckDataMerge()`, `check_srt_list()` to `CheckDataList` and `check_data_type()` to `CheckDataType`.

# scop 0.2.2

* **feat**:
  * Replace all `BiocParallel::bplapply()` with `thisutils::parallelize_fun()`.

* **fix**:
  * Fix bugs in `RunSingleR()`.

# scop 0.2.0

* **feat**:
  * Added `remove_python()` function for easy remove *Python* packages.

* **fix**:
  * Corrected an issue in `py_to_r2()` function (intrinsic function), which ensures that Python-dependent functions like `RunPAGA()` and `RunSCVELO()` function run correctly.

# scop 0.1.9

* **feat**:
  * Update `CellScoring()` and `AddModuleScore2()` functions. Now, new parameters `cores` and `verbose` have been added.
  * `AddModuleScore2()` function no longer uses the `BiocParallel::bpparam()` function to enable parallelization, but `thisutils::parallelize_fun`, and the `cores` parameter is used to control the number of cores in `thisutils::parallelize_fun`.

# scop 0.1.5

* **fix**:
  * Fix error for `RunPalantir()` function, see [#23](https://github.com/mengxu98/scop/issues/23).

# scop 0.1.4

* **feat**:
  * Update `.onAttach()`, now `.onAttach()` will print more information about *conda* and *Python*.
  * Update `PrepareEnv()` function for easy add or update a conda environments and install *Python* packages.
  * Added `ListEnv()` and `RemoveEnv()` functions for easy management of *conda* environment and *Python* packages.

# scop 0.1.3

* **feat**:
  * Added `TACSPlot()` function for creating FACS-like plots. Please refer to [Kernfeld et al. paper](https://doi.org/10.1016/j.immuni.2018.04.015) and [Github](https://github.com/maehrlab/thymusatlastools2/blob/f8b51ad684d56b2eeda780787eb9ad4ff3003eef/R/data_handling_seurat.R#L271) for specific information.

# scop 0.0.9

* **fix**:
  * Fix a bug for `PrepareEnv()` function.
  * The default *Python* version is now set to `3.10-1`.

# scop 0.0.6

* **feat**:
  * Add `GetAssayData5()` function, a reimplementation of `GetAssayData()`, for compatibility with Seurat v5 `Assay` objects.
  * Updated `GetFeaturesData()` and `AddFeaturesData()` function to support retrieving and adding feature metadata.

# scop 0.0.5

* **data**:
  * Updated the `pancreas_sub` and `panc8_sub` test datasets to the Seurat v5 object format.

# scop 0.0.1

* **Initial version**
