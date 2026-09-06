# Changelog

## scop 0.9.1

- **fix**: Coordinate v3 accepts the empty native orientation slot
  introduced by newer SeuratObject versions. List integration resolves
  explicit image maps before merging duplicate image names. Bundled
  metadata-only result fixtures are migrated only after pinned-input
  coordinate-equivalence and payload-preservation checks; arbitrary user
  v2 results still require rerunning.

- **fix**: Spatial coordinate contract v3 uses persistent image axis
  provenance instead of guessing from mutable metadata, parses numeric
  factor coordinates correctly, and records the resolved coordinate
  columns.
  [`SetSpatialImageAxes()`](https://mengxu98.github.io/scop/reference/SetSpatialImageAxes.md)
  preserves custom VisiumV2 conventions through subsetting and renaming;
  the bundled Visium object is migrated without changing its
  measurements. Coordinate-dependent results from contract v2 must be
  rerun.

- **fix**: Spatial point and pie plots default to the same FOV/metadata
  y orientation as network and boundary plots (`flip.y = FALSE`).
  SpatialExperiment exports default to raw coordinates and retain
  transforms for reversible display-coordinate round trips.

- **fix**: Spatial integration and SpatialEcoTyper multi discovery
  resolve each sample’s image independently and reject incomplete image
  coverage. SpatialEcoTyper single discovery requires a single selected
  image and stores its raw coordinate source. SpotSweeper local outlier
  results with verified sample-prefixed IDs are restored to input order
  without accepting changed coordinates or expression.

- **fix**: Gradient result storage and VariableFeatures updates are
  independent; explicit empty updates clear VariableFeatures and
  completion messages describe actual actions.

- **fix**: Neighborhood spatial counts respect saved analysis scope and
  condition, distinguish zero from unevaluated cells, reject unsupported
  plot arguments, and support valid empty observed results.

- **feat**:

  - [`SpatialNeighborhoodProfile()`](https://mengxu98.github.io/scop/reference/SpatialNeighborhoodProfile.md)
    returns observed per-cell neighbor counts and fractions at multiple
    raw-coordinate distances, retaining the full selected tissue context
    for target cells without modifying the Seurat object or storing an
    edge table.
  - [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md)
    sample-level methods gain `min.cells.sample` (default `10`) to drop
    pseudobulk samples whose cell count in the current `group.by` level
    is below a threshold, which can reverse DE direction relative to
    cell-level tests. Set `min.cells.sample = 1` to keep the previous
    unfiltered behavior.
  - Internal preprocessing and integration workflows now use scop’s
    Seurat-compatible entry points, enabling validated native
    acceleration for normalization, variable-feature selection, scaling,
    PCA, and neighbor search while retaining transparent Seurat fallback
    for unsupported paths.
  - The native `SCTransform(vst.flavor = "v2")` path now supports
    validated numeric and factor `vars.to.regress` designs and
    `latent.data`; unsupported regression designs continue to delegate
    to Seurat.
  - The native PCA auto path accepts smaller stable trailing eigengaps
    after numerical and neighborhood-parity validation, allowing more
    datasets to use the accelerated implementation.
  - Added
    [`CellScoringPlot()`](https://mengxu98.github.io/scop/reference/CellScoringPlot.md)
    for paired group/score UMAPs, score distributions, threshold
    proportions, and complete multi-program layouts from
    [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md)
    results.
  - CellRank + Palantir fate workflow now supports MAGIC layers,
    explicit Schur sizes and terminal-state validation, normalized
    fate/driver/transition payloads, CellRank GAM trend modules, module
    enrichment, and result-only
    [`CellRankPlot()`](https://mengxu98.github.io/scop/reference/CellRankPlot.md)
    views.
  - Added
    [`RunCOMMOT()`](https://mengxu98.github.io/scop/reference/RunCOMMOT.md)
    and
    [`COMMOTPlot()`](https://mengxu98.github.io/scop/reference/COMMOTPlot.md)
    for official COMMOT spatial communication through an isolated Python
    subprocess with raw-coordinate provenance and optional external H5AD
    storage.
  - Added
    [`RunSpaTalk()`](https://mengxu98.github.io/scop/reference/RunSpaTalk.md)
    and
    [`SpaTalkPlot()`](https://mengxu98.github.io/scop/reference/SpaTalkPlot.md)
    for official SpaTalk single-cell and spot-reference spatial
    communication workflows with compact unified CCC storage.
  - Added `CCCNetworkPlot(plot_type = "spatial")`, a stored-result slice
    overlay for SpatialCellChat, SpaTalk, and COMMOT with
    abundance-scaled group nodes, sender-coloured directed edges, and
    proportion-aware composition pies.
  - [`RunCNV()`](https://mengxu98.github.io/scop/reference/RunCNV.md)
    backend coverage: `copykat`, `fastCNV`, `SCEVAN`, `infercnv`, and
    `numbat` are now verified against their original pipelines on real
    bundled datasets (breast-cancer scRNA-seq, infercnv example, numbat
    ATC2 example), with exact prediction agreement for
    copykat/SCEVAN/numbat and Spearman rho = 1.0 for fastCNV/infercnv.
  - The MERINGUE C++ permutation kernel now traverses the weight matrix
    sparsely (row-major edge list, bit-identical accumulation):
    `nperm = 200` Moran p-values at 300 spots drop from ~442s to ~14s
    (~30x).
  - [`RunMERINGUE()`](https://mengxu98.github.io/scop/reference/RunMERINGUE.md)
    gains a `backend = c("cpp", "r")` option; the C++ kernels replicate
    `MERINGUE::moranTest()`/`moranPermutationTest()` to machine
    precision (including the `sample.kind` RNG stream). `backend = "r"`
    retains the original calls and honors `moran_params`.
  - Real-data consistency tests verify the wrappers against the original
    pipelines (CellChat, SpatialCellChat, LIANA, BayesSpace, MERINGUE,
    SPOTlight, SpaNorm, SpatialDWLS, RCTD, BANKSY, CARD, SpotSweeper,
    PRECAST).
  - [`DorotheaPlot()`](https://mengxu98.github.io/scop/reference/DorotheaPlot.md)
    now follows the shared scop plot style: horizontal
    `"bar"`/`"lollipop"` comparison plots, `"volcano"`, `"heatmap"` via
    [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md),
    `"dim"` (activity vs expression on embeddings), `"stat"`
    (`FeatureStatPlot` distributions), and `"targets"` (regulon-target
    volcano).
    [`RunDorothea()`](https://mengxu98.github.io/scop/reference/RunDorothea.md)
    examples show
    [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
    /
    [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)
    on the TF-activity assay.
  - [`RunIntegrationBenchmark()`](https://mengxu98.github.io/scop/reference/RunIntegrationBenchmark.md)
    compares integration methods with scIB-style iLISI/cLISI, ASW, graph
    connectivity, and ARI/NMI scores and returns a `Seurat` object.
    Tables live in `srt@tools$IntegrationBenchmark` and are printed with
    [`thisplot::print_colored_table()`](https://mengxu98.github.io/thisplot/reference/print_colored_table.html).
    [`IntegrationBenchmarkPlot()`](https://mengxu98.github.io/scop/reference/IntegrationBenchmarkPlot.md)
    draws `box` (per-cell LISI), heatmap, scatter, and UMAP views.
    `LISIPlot()` is removed; use
    `IntegrationBenchmarkPlot(plot_type = "box")`.
  - [`SCENICPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlot.md)
    /
    [`SCENICPlusPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlusPlot.md):
    `heatmap_dotplot`, `eregulon_dim`, `coverage`, `network`,
    `network_graph`, `egrn`, `overlap`. `"network"` draws one hub per
    TF. Network plots have no title. Networks use `"Chinese"` when
    `palette = "RdYlBu"`. Multi-TF legends show the top RSS cell type.

- **fixed**:

  - C-SIDE completion messages respect result storage and document
    compatibility metadata as run-level significant record counts.
  - `presto` is now a runtime-optional GitHub backend instead of a
    suggested check dependency. Marker fast paths use an existing
    installation without installing it, while
    `RunCellChat(do.fast = TRUE)` can resolve the backend from
    `immunogenomics/presto` before execution.
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md)
    now disables PETSc explicitly when retrying fate probabilities with
    SciPy’s direct solver, preventing CellRank 2.0.7 from silently
    reverting the retry to GMRES and rejecting harmless floating-point
    negative probabilities
    ([\#380](https://github.com/mengxu98/scop/issues/380)).
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md)
    defaults to the dependency-free `brandts` Schur method again, so the
    standard `scop_env` no longer fails because the optional `petsc4py`
    and `slepc4py` packages are absent. Users can still request
    `schur_method = "krylov"` when those packages are installed
    ([\#380](https://github.com/mengxu98/scop/issues/380)).
  - CellRank stored-result plotting now follows SCOP’s `theme_scop`,
    `Chinese`, and `Spectral` visual contracts; draws transition
    projections directly from the saved sparse transition matrix;
    supports random walks without requiring trend results; facets trend
    modules; and adds module-enrichment dot plots. CellRank enrichment
    selects the database gene-ID column that overlaps the actual driver
    background, including three-column GO tables with both Entrez and
    symbol identifiers.
  - Linux shared-library preloading keeps a package-local
    `load_shared_library()` seam, so dependency-order retries can be
    tested without attempting to mock the base
    [`dyn.load()`](https://rdrr.io/r/base/dynload.html) function as a
    SCOP namespace binding.
  - Coordinate-contract validation now covers unified CCC plots for
    SpatialCellChat, COMMOT, and SpaTalk without rejecting
    coordinate-independent SpatialDWLS proportions. Explicit
    `coordinate_space = "legacy_display"` analysis preserves its
    historical lowres-to-hires fallback; current display plotting
    remains strict about the selected scale.
  - Spatial coordinate contract v2 now reads full-resolution
    VisiumV1/VisiumV2 coordinates, restores the `imagerow`/`imagecol`
    order that Seurat’s standard VisiumV2 loader stores positionally as
    centroid `x`/`y`, preserves explicit generic FOV `x`/`y`, normalizes
    only explicitly named legacy image row/column coordinates, and
    applies exactly one selected `lowres` or `hires` display transform.
    [`SpatialCoordinates()`](https://mengxu98.github.io/scop/reference/SpatialCoordinates.md),
    HE-overlay point/pie/network plots, and their public wrappers gain
    an appended `image.scale` argument; variadic wrappers keep it
    named-only after `...`. Raw SpatialCellChat analysis no longer
    requires an available display scale, and coordinate-independent
    producers no longer claim a coordinate-contract version. Users no
    longer need to modify Seurat `@scale.factors` slots.
    Coordinate-dependent stored results created under an older or
    unknown contract must be rerun.
  - [`RunSlingshot()`](https://mengxu98.github.io/scop/reference/RunSlingshot.md)
    / `CellDimPlot(lineages=)` no longer fail with
    `scale$guide == "none"` on ggplot2 \>= 3.5 when `ggnewscale` is
    older than 0.5.0. Colour/fill scales now keep a character `guide`,
    legend overrides go through `guides()`, `ggnewscale` is pinned to
    \>= 0.5.2, and a failed lineage plot no longer discards the
    Slingshot result.
  - [`RunGSEA()`](https://mengxu98.github.io/scop/reference/RunGSEA.md)
    /
    [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md)
    no longer fail with “replacement has 0 rows, data has N” when the
    cached database carries an empty version. The KEGG REST API dropped
    the `Release` line from `/info/<org>`, so
    [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md)
    now resolves the KEGG version from the release-notes page when the
    info endpoint does not report it (retrieval-date fallback; HPO
    hardened the same way), and the enrichment writers tolerate missing
    versions.
  - [`RunSpaTalk()`](https://mengxu98.github.io/scop/reference/RunSpaTalk.md)
    now reuses stored
    [`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md)
    weights through SpaTalk method 2, records the deconvolution source,
    and accepts a working installed SpaTalk package even when its
    DESCRIPTION lacks GitHub remote metadata.
  - [`RunLIANA()`](https://mengxu98.github.io/scop/reference/RunLIANA.md)
    now preserves the field set returned by LIANA’s official consensus
    mode: rank consensus no longer gains all-`NA` aggregate/significance
    placeholders, aggregate consensus no longer gains all-`NA`
    magnitude/specificity placeholders, and deterministic
    ligand-receptor display identifiers are completed during SCOP/LIANA
    export.
  - [`RunCNV()`](https://mengxu98.github.io/scop/reference/RunCNV.md)
    attaches backend packages to the search path during execution so
    lazy data objects (copykat `full.anno`, fastCNV `geneMetadata`,
    numbat `gtf_hg38`) are resolvable from backend namespace code;
    `copykat` is installed from its GitHub source
    (`navinlabcode/copykat`) because it was removed from CRAN.
  - `RunCNV(numbat)` prefers the plain `seg` column over `seg_label`
    when reading `joint_post` output, so segment identifiers are not
    polluted with the `(loh)` suffix.
  - [`RunCNV()`](https://mengxu98.github.io/scop/reference/RunCNV.md)
    numbat extraction is vectorized (direct index assignment with
    duplicate averaging instead of row-by-row `aggregate`/`strsplit`,
    plus [`split()`](https://rdrr.io/r/base/split.html)-partitioned bin
    statistics): a 300-seg x 2000-cell `joint_post` extraction drops
    from ~16s to ~1.3s (~12x) with bit-identical results.
  - `visium_mouse_brain_slices_sub` rebuilt with the current
    `SeuratObject` so bundled `VisiumV2` images can be subset again.
  - [`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md)
    extracts per-spot `all_weights` from `spacexr` \>= 2.2
    `doublet_mode = "multi"` results.
  - [`RunSpotSweeper()`](https://mengxu98.github.io/scop/reference/RunSpotSweeper.md)
    detects upstream
    [`SpotSweeper::localVariance`](https://rdrr.io/pkg/SpotSweeper/man/localVariance.html)
    breakage on `SpatialExperiment` colData and degrades gracefully.
  - `RunSpatialIntegration(PRECAST)` writes the raw coordinates expected
    by `CreatePRECASTObject`.
  - [`RunLIANA()`](https://mengxu98.github.io/scop/reference/RunLIANA.md)
    recognizes `logfc_comb` (and other liana method score columns) as
    the interaction score.
  - Exported `theme_scop` / `theme_spatial` are listed in the pkgdown
    reference index; the `theme_scop` alias test compares function
    equality rather than identity so installed-package checks pass.
  - `DorotheaPlot(plot_type = "volcano")` plots every tested TF, caps
    Wilcoxon underflow on the y-axis, colors points by direction, and
    labels significant TFs. `"targets"` labels only significant genes so
    leader lines no longer pile at x = 0.

- **removed**:

  - Removed wrappers whose upstream dependencies are unavailable or
    uncompilable: `RunSpatialQM()`, the `RunSemla*()` wrappers, the
    `BASS`/`SpatialMNN` branches of
    [`RunSpatialIntegration()`](https://mengxu98.github.io/scop/reference/RunSpatialIntegration.md)
    (PRECAST unaffected), and the SPATA2 `backend = "r"` path of
    [`RunSpatialGradientFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialGradientFeatures.md)
    (C++ backend only; `srt_to_spata2()`/`spata2_to_srt()` removed as
    well).
  - Removed the legacy spatial registry and result-discovery APIs
    (`GetSpatialResult()`, `ListSpatialMethods()`,
    `SpatialBackendStatus()`, and `SpatialResultInfo()`). Use each
    producer’s documented `srt@tools` entry and, where provided, its
    dedicated result accessor or plotting function.
  - The `CaSpER` branch of
    [`RunCNV()`](https://mengxu98.github.io/scop/reference/RunCNV.md) is
    removed: `runCaSpER()` segfaults deterministically on R \>= 4.5 (HMM
    `manualSegment`, even on 50 cells), upstream is unmaintained since
    2019, and its Bioconductor dependency set is archived.
    [`RunCNV()`](https://mengxu98.github.io/scop/reference/RunCNV.md)
    now supports `copykat`, `fastCNV`, `scevan`, `infercnv`, and
    `numbat`.

- **changed**:

  - `RunBenchmark()` / `BenchmarkPlot()` are replaced by
    [`RunSpatialBenchmark()`](https://mengxu98.github.io/scop/reference/RunSpatialBenchmark.md)
    /
    [`SpatialBenchmarkPlot()`](https://mengxu98.github.io/scop/reference/SpatialBenchmarkPlot.md).
    The old names are removed.
  - [`RunCell2fate()`](https://mengxu98.github.io/scop/reference/RunCell2fate.md)
    examples call
    [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md)
    and plot time, uncertainty, module activation, and module state.
    [`env_requirements()`](https://mengxu98.github.io/scop/reference/env_requirements.md)
    also shows the standalone Cell2fate Python 3.9 stack.
  - `ListLIANAResources()` is replaced by
    [`ListCCCDB()`](https://mengxu98.github.io/scop/reference/ListCCCDB.md),
    which enumerates ligand-receptor databases and prior models across
    all CCC backends with a unified `db`/`species`/`status` schema.
  - New
    [`PrepareCCCDB()`](https://mengxu98.github.io/scop/reference/PrepareCCCDB.md)
    prepares the CellTalk and CellChat databases
    (`TERM2GENE`/`TERM2NAME`, `R.cache`-backed); the CellTalk/CellChat
    branches were removed from
    [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md).
  - Database listing outputs are unified
    ([`ListDB()`](https://mengxu98.github.io/scop/reference/ListDB.md)
    and
    [`ListCCCDB()`](https://mengxu98.github.io/scop/reference/ListCCCDB.md)
    return plain data frames; colored console rendering is available via
    [`thisplot::print_colored_table()`](https://mengxu98.github.io/thisplot/reference/print_colored_table.html)).
  - Python-converted object interop:
    [`PAGAPlot()`](https://mengxu98.github.io/scop/reference/PAGAPlot.md),
    [`VelocityPlot()`](https://mengxu98.github.io/scop/reference/VelocityPlot.md)/`CellDimPlot(velocity=)`,
    `RunPAGA(use_rna_velocity)`,
    [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md),
    DE consumers
    ([`VolcanoPlot()`](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
    [`DEtestPlot()`](https://mengxu98.github.io/scop/reference/DEtestPlot.md),
    [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
    [`RunGSEA()`](https://mengxu98.github.io/scop/reference/RunGSEA.md),
    …), and
    [`RunScissor()`](https://mengxu98.github.io/scop/reference/RunScissor.md)
    now fall back to scanpy/scVelo result slots (`paga`,
    `rank_genes_groups`, `velocity_umap`, `dpt_pseudotime`,
    `connectivities`) stored by
    [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md).

## scop 0.9.0

- **feat**:
  - New optional-backend wrappers:
    [`RunCHOIR()`](https://mengxu98.github.io/scop/reference/RunCHOIR.md)
    for significance-based hierarchical clustering and
    [`RunCell2fate()`](https://mengxu98.github.io/scop/reference/RunCell2fate.md)
    for RNA-velocity modeling, both with pinned backend provenance,
    strict result validation, metadata writeback, and resumable
    artifacts.
  - Added
    [`RunCell2location()`](https://mengxu98.github.io/scop/reference/RunCell2location.md)
    for official Python cell2location reference-signature learning and
    spatial abundance mapping, with manifest-validated resume, signature
    reuse, and
    [`Cell2locationPlot()`](https://mengxu98.github.io/scop/reference/Cell2locationPlot.md)
    visualizations.
  - Added
    [`RunCNV()`](https://mengxu98.github.io/scop/reference/RunCNV.md) as
    a unified CNV workflow with optional `copykat`, `infercnv`,
    `SCEVAN`, and `fastCNV` backends, standardized `srt@tools[["CNV"]]`
    storage, and
    [`CNVPlot()`](https://mengxu98.github.io/scop/reference/CNVPlot.md)
    visualizations; added
    [`RunESTIMATE()`](https://mengxu98.github.io/scop/reference/RunESTIMATE.md)
    for stromal/immune/tumor-purity scoring with
    [`EstimateScorePlot()`](https://mengxu98.github.io/scop/reference/EstimateScorePlot.md)
    /
    [`EstimateGenePlot()`](https://mengxu98.github.io/scop/reference/EstimateGenePlot.md).
  - GRN analysis extended:
    [`RunSCENICPlus()`](https://mengxu98.github.io/scop/reference/RunSCENICPlus.md)
    for SCENIC+ multi-omics workflows, plus standalone
    [`RunGRNBoost2()`](https://mengxu98.github.io/scop/reference/RunGRNBoost2.md),
    [`RunGENIE3()`](https://mengxu98.github.io/scop/reference/RunGENIE3.md),
    and
    [`RunCisTarget()`](https://mengxu98.github.io/scop/reference/RunCisTarget.md)
    modules sharing the SCENIC regulon infrastructure.
  - Added
    [`RunSpatialCellChat()`](https://mengxu98.github.io/scop/reference/RunSpatialCellChat.md)
    as a SpatialCellChat v3 producer for cell-, spot-, and
    composition-level communication, with raw-to-micron coordinate
    contracts and schema-v1 provenance;
    [`GetCCCObject()`](https://mengxu98.github.io/scop/reference/GetCCCObject.md)
    and
    [`SpatialCellChatPlot()`](https://mengxu98.github.io/scop/reference/SpatialCellChatPlot.md)
    provide stored-result access.
  - Spatial analysis framework: added spatial contract schema v1
    ([`SpatialCoordinates()`](https://mengxu98.github.io/scop/reference/SpatialCoordinates.md),
    `GetSpatialResult()`,
    [`GetSpatialGraph()`](https://mengxu98.github.io/scop/reference/GetSpatialGraph.md)),
    a shared spatial method registry (`ListSpatialMethods()`,
    `SpatialBackendStatus()`, `SpatialResultInfo()`), strict
    raw-coordinate defaults, and sparse `BiocNeighbors`-backed graph
    kernels for SmoothClust, spatial-variable-feature, and
    spatial-neighborhood paths.
  - Added
    [`RunSpatialBenchmark()`](https://mengxu98.github.io/scop/reference/RunSpatialBenchmark.md)
    and
    [`SpatialBenchmarkPlot()`](https://mengxu98.github.io/scop/reference/SpatialBenchmarkPlot.md)
    for failure-tolerant comparison of stable spatial-domain clustering
    methods (ARI/NMI/purity, runtime, peak memory); added
    [`ClusterTreePlot()`](https://mengxu98.github.io/scop/reference/ClusterTreePlot.md)
    for multi-resolution clustering trees.
  - Visualization improvements:
    [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)/`ExpressionStatPlot()`
    gain `auto_comparison`;
    [`GSVAPlot()`](https://mengxu98.github.io/scop/reference/GSVAPlot.md)
    gains a `mode = "diff"` two-group pathway test;
    [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md),
    [`DynamicHeatmap()`](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md),
    and
    [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)
    compact long multi-term annotations;
    [`CellCorHeatmap()`](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md)
    gains `legend.position`.
  - Performance: C++ acceleration for the core Seurat flow —
    [`NormalizeData()`](https://mengxu98.github.io/scop/reference/NormalizeData.md),
    [`FindVariableFeatures()`](https://mengxu98.github.io/scop/reference/FindVariableFeatures.md),
    [`ScaleData()`](https://mengxu98.github.io/scop/reference/ScaleData.md),
    PCA and neighbor search — and for marker detection
    ([`FindMarkers()`](https://mengxu98.github.io/scop/reference/FindMarkers.md),
    [`FindAllMarkers()`](https://mengxu98.github.io/scop/reference/FindAllMarkers.md),
    `FindConservedMarkers()`); native kernels for manifold learning
    (UMAP, PaCMAP, TriMap, PHATE, DM), trajectory inference (PAGA,
    SCVELO, VECTOR), Z-score/PLAGE/ssGSEA/AUCell gene-set scoring,
    spatial-variable-feature statistics, cross-KNN and connectivity
    graphs, and integration embeddings.
  - [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md)
    now builds gene sets via
    [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md)
    by default with species-aware conversion;
    [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md)
    accepts local single-file database sources (MSigDB, CORUM, JASPAR,
    TRRUST, CellChat, and more);
    [`RunSCENIC()`](https://mengxu98.github.io/scop/reference/RunSCENIC.md)
    gained a `genome` parameter for hg38/hg19 reference selection.
  - [`RunSpatialNeighborhood()`](https://mengxu98.github.io/scop/reference/RunSpatialNeighborhood.md)
    defaults to observed-neighborhood summaries; the standard spatial
    workflow records per-stage requested/completed/skipped/failed state
    and reports partial completion.
- **changed**:
  - Renamed
    [`standard_scop()`](https://mengxu98.github.io/scop/reference/deprecated-workflows.md)
    to
    [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md)
    and
    [`integration_scop()`](https://mengxu98.github.io/scop/reference/deprecated-workflows.md)
    to
    [`RunIntegration()`](https://mengxu98.github.io/scop/reference/RunIntegration.md)
    (multi-method support via `integration_methods`); renamed
    [`db_scDblFinder()`](https://mengxu98.github.io/scop/reference/deprecated-doublet-callers.md),
    [`db_scds()`](https://mengxu98.github.io/scop/reference/deprecated-doublet-callers.md),
    [`db_Scrublet()`](https://mengxu98.github.io/scop/reference/deprecated-doublet-callers.md),
    and
    [`db_DoubletDetection()`](https://mengxu98.github.io/scop/reference/deprecated-doublet-callers.md)
    to
    [`RunscDblFinder()`](https://mengxu98.github.io/scop/reference/RunscDblFinder.md),
    [`Runscds()`](https://mengxu98.github.io/scop/reference/Runscds.md),
    [`RunScrublet()`](https://mengxu98.github.io/scop/reference/RunScrublet.md),
    and
    [`RunDoubletDetection()`](https://mengxu98.github.io/scop/reference/RunDoubletDetection.md).
    Old entry points emit deprecation warnings and will be removed in
    version 1.0.0.
  - Spatial analysis, plotting, and framework conversion now require
    explicit `image` selection for multi-image objects;
    distance-sensitive producers default to raw acquisition coordinates
    (`coordinate_space = "legacy_display"` remains available).
- **data**:
  - Added `visium_human_pancreas_sub`, a Visium human pancreas spatial
    example dataset
    ([GSE254829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254829))
    with spatial image data and CODA-derived annotations, for spatial
    analysis.
  - Added a Xenium package-data import script for preparing curated
    Xenium example resources outside the installed package payload.

## scop 0.8.9

- Cell-cell communication:
  - [`RunLIANA()`](https://mengxu98.github.io/scop/reference/RunLIANA.md)
    now stores LIANA’s official multi-method consensus from
    `rank_aggregate()` or `liana_aggregate()` alongside the historical
    method-stacked and scop-postprocessed tables. Human and mouse
    defaults use the `Consensus` and `MouseConsensus` ligand-receptor
    resources, respectively, and multiple resources are aggregated
    independently.
  - Added `ListLIANAResources()`, `ListCCCMethods()`, `CCCResultInfo()`,
    and `GetCCCResult()` for resource discovery and Seurat-centered CCC
    result access.
  - [`RunCCC()`](https://mengxu98.github.io/scop/reference/RunCCC.md)
    now registers all seven SCOP CCC producers while retaining CellChat,
    CellPhoneDB, and LIANA as the default core run. Unified CCC plotting
    separates backends by default and provides explicit support-count
    and within-method percentile-rank views.
  - Hardened unified CCC plotting and access semantics: context filters
    are applied before cross-method support/rank combinations, CellChat
    cached rows honor the requested significance threshold, LIANA
    resources remain separate during export, and unified retrieval
    preserves backend labels. `RunCCC(skip_failed = TRUE)` now records
    preflight failures without aborting otherwise runnable methods.
- **fix**:
  - [`RunGSVA()`](https://mengxu98.github.io/scop/reference/RunGSVA.md):
    Removed redundant dense
    [`as.matrix()`](https://rdrr.io/r/base/matrix.html) conversion in
    single-cell mode for `backend = "r"`; sparse expression matrices now
    stay sparse through row filtering, reducing peak memory and avoiding
    unnecessary dense materialization.
  - [`RunCellQC()`](https://mengxu98.github.io/scop/reference/RunCellQC.md):
    Replaced
    [`Seurat::SplitObject()`](https://satijalab.org/seurat/reference/SplitObject.html)
    with lazy cell-name-based subsetting, avoiding full-object
    duplication per split group. Added caching of `GetAssay()` gene
    names,
    [`GetAssayData5()`](https://mengxu98.github.io/scop/reference/GetAssayData5.md)
    counts, and per-cell QC metrics (`nCount` / `nFeature`) outside the
    species-check loop to eliminate repeated data extraction.
  - [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md):
    Pre-computed cell index mapping
    (`split(names(cell_group), cell_group)`) for paired-marker tests,
    replacing per-pair [`which()`](https://rdrr.io/r/base/which.html)
    calls with O(1) list lookups.
  - [`AnnotateFeatures()`](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md):
    Replaced per-detail `sapply(… "[")` with type-stable
    `vapply(…, character(1))` to avoid implicit list-to-character
    coercion.
  - `run_scomm()`: Deferred dense conversion for reference and query
    subsetting by removing premature
    [`as.matrix()`](https://rdrr.io/r/base/matrix.html) calls; sparse
    matrices are subset first, then densified only when required.
  - [`RunUMAP2()`](https://mengxu98.github.io/scop/reference/RunUMAP2.md):
    Replaced `isSymmetric(as.matrix(graph))` with sparse-aware
    `Matrix::isSymmetric(graph)` in symmetry checks and graph subset
    sampling, avoiding dense materialization of large neighbor graphs.
  - [`RunUMAP2()`](https://mengxu98.github.io/scop/reference/RunUMAP2.md):
    Replaced `apply(as.matrix(graph), 2, order)` in the uwot-predict
    path with `thisutils::run_sparse_topk_stored(by = "col")`, avoiding
    full dense conversion and column-wise R-level
    [`apply()`](https://rdrr.io/r/base/apply.html) while retaining
    graph-edge semantics.
  - [`RunDimsReduction()`](https://mengxu98.github.io/scop/reference/RunDimsReduction.md)
    (PCA centering): Uses
    `SeuratObject::LayerData(…, features = features)` to read only HVF
    rows from `scale.data` instead of loading the full matrix followed
    by manual subsetting.
  - [`RunMDS()`](https://mengxu98.github.io/scop/reference/RunMDS.md):
    Removed [`as.matrix()`](https://rdrr.io/r/base/matrix.html) before
    `Matrix::t()`;
    [`proxyC::dist`](https://koheiw.github.io/proxyC/reference/simil.html)
    now receives the sparse matrix directly, avoiding dense conversion
    of large count matrices.
  - `integration.R` (fastMNN): Removed three
    [`as.matrix()`](https://rdrr.io/r/base/matrix.html) calls passed to
    [`batchelor::fastMNN()`](https://rdrr.io/pkg/batchelor/man/fastMNN.html),
    which accepts sparse matrices directly via `SingleCellExperiment`.
  - [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md):
    Replaced full `scale.data` matrix load (via
    [`GetAssayData5()`](https://mengxu98.github.io/scop/reference/GetAssayData5.md))
    with direct layer access
    ([`SeuratObject::GetAssayData`](https://satijalab.github.io/seurat-object/reference/AssayData.html)
    for Assay5; `@scale.data` for Assay) to retrieve only rownames when
    checking whether HVFs have been scaled, substantially reducing peak
    memory during the ScaleData decision step.
  - [`GetAssayData5.Assay()`](https://mengxu98.github.io/scop/reference/GetAssayData5.md):
    Fixed parameter naming to use positional matching for the slot/layer
    argument, so
    [`GetAssayData5()`](https://mengxu98.github.io/scop/reference/GetAssayData5.md)
    correctly retrieves `counts`, `data`, and `scale.data` layers in
    both Seurat v4 (`slot`) and v5 (`layer`). Previously the named
    `layer` argument was silently dropped by SeuratObject v4, always
    returning the `data` slot regardless of the requested layer.
  - [`CSS_integrate()`](https://mengxu98.github.io/scop/reference/CSS_integrate.md):
    Added `Assay5` guard before
    [`SeuratObject::JoinLayers()`](https://satijalab.github.io/seurat-object/reference/SplitLayers.html)
    to avoid errors with Seurat v4 `Assay` objects, which do not support
    layered storage.
  - [`RunDimsReduction()`](https://mengxu98.github.io/scop/reference/RunDimsReduction.md):
    Added Seurat v4 fallback for PCA centering —
    `SeuratObject::LayerData(…, features = …)` is Assay5-only; v4
    `Assay` objects now use
    [`GetAssayData5()`](https://mengxu98.github.io/scop/reference/GetAssayData5.md)
    with manual feature subsetting.
  - [`RunIntegration()`](https://mengxu98.github.io/scop/reference/RunIntegration.md):
    Added early detection for v5-only integration methods (`CCA`,
    `RPCA`, `fastMNN5`, `Harmony5`, `scVI5`) to provide a clear,
    actionable error message when used with Seurat v4 `Assay` objects,
    rather than failing deep in the call stack.
- **feat**:
  - Added
    [`loom_to_srt()`](https://mengxu98.github.io/scop/reference/loom_to_srt.md)
    for pure-R loom-to-Seurat conversion via `rhdf5`, preserving
    velocity `spliced` and `unspliced` layers as assays without
    initializing Python, and added Python-backed
    [`loom_to_adata()`](https://mengxu98.github.io/scop/reference/loom_to_adata.md)
    for users who need AnnData output.
  - Added `RunBulk()` as a unified bulk-strategy entrypoint with
    method-vector selection for bulk DE, deconvolution, and
    cell-type-specific DE workflows.
  - Added method-specific bulk runners for `de_limma_voom`,
    `de_edgeR_qlf`, `de_DESeq2`, `de_dream`, `deconv_MuSiC`,
    `deconv_BisqueRNA`, `deconv_BayesPrism`, and `csde_TOAST`.
  - Standardized bulk results under `Bulk$results$de`,
    `Bulk$results$deconv`, and `Bulk$results$csde`, keeping DE outputs
    compatible with existing
    [`DEtestPlot()`](https://mengxu98.github.io/scop/reference/DEtestPlot.md),
    [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
    [`RunGSEA()`](https://mengxu98.github.io/scop/reference/RunGSEA.md),
    [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md),
    and
    [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md)
    data flows.
  - `RunBulk(run_enrichment = TRUE, run_gsea = TRUE)` now filters bulk
    DE rows with pathway thresholds and mirrors successful pathway
    results to `Enrichment_Bulk_wilcox` and `GSEA_Bulk_wilcox`, so
    [`EnrichmentPlot()`](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md)
    and
    [`GSEAPlot()`](https://mengxu98.github.io/scop/reference/GSEAPlot.md)
    can read them through the standard `srt@tools` contract.
  - Deconvolution and CSDE bundles now record their computational
    `engine` in `details`; the current deconvolution runners use
    explicit `backend = "internal"` SCOP profile fitting and
    `csde_TOAST` uses explicit `backend = "limma_interaction"`, so
    package backends can be wired in later without pretending to call
    external packages.
  - Added explicit `Remotes` entries for `mengxu98/thisplot` and
    `mengxu98/thisutils` so source installs can resolve the minimum
    imported versions required by SCOP. Optional bulk engines such as
    `DESeq2` and `variancePartition` are installed on demand through
    `check_r()` and called lazily through `get_namespace_fun()`.
  - Added `method_args` to expose method-specific tuning parameters
    without expanding the public API surface.

## scop 0.8.8

- **deps**:
  - Updated the minimum dependency versions to `thisplot (>= 0.3.8)` and
    `thisutils (>= 0.4.5)`, and removed local copies of helpers now
    provided upstream
    ([`thisplot::annotate_quadrants()`](https://mengxu98.github.io/thisplot/reference/annotate_quadrants.html),
    [`thisplot::clip_symmetric_range()`](https://mengxu98.github.io/thisplot/reference/clip_symmetric_range.html),
    [`thisutils::collapse_sparse_rows()`](https://mengxu98.github.io/thisutils/reference/collapse_sparse_rows.html)).
- **docs**:
  - Updated pkgdown reference grouping for the cell-cycle workflow.
- **fix**:
  - Added optional wrappers for
    [`RunDorothea()`](https://mengxu98.github.io/scop/reference/RunDorothea.md),
    [`RunBayesSpace()`](https://mengxu98.github.io/scop/reference/RunBayesSpace.md),
    and experimental
    [`RunscTenifoldKnk()`](https://mengxu98.github.io/scop/reference/RunScTenifoldKnk.md).
    [`RunscTenifoldKnk()`](https://mengxu98.github.io/scop/reference/RunScTenifoldKnk.md)
    keeps the upstream `scTenifoldNet` workflow but fixes the QC
    gene-filter assignment in the local path, uses an equivalent
    covariance/downdate path with direct sparse matrix construction,
    selection-based quantile thresholding, and controlled per-gene
    eigensolver parallelism for large `pcNet()` network construction,
    and uses C++ helpers for tensor decomposition, manifold matrix
    construction, directionality, and differential-regulation distance
    calculations. The C++ tensor-decomposition path now computes MTTKRP
    updates directly instead of materializing four dense unfolding
    matrices.
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    Added explicit support for `mamba` and `micromamba` executables in
    addition to `conda`, including command-name resolution, automatic
    package-managed micromamba download when `conda = "micromamba"` is
    requested and the command is not on `PATH`, manager-aware logging,
    micromamba env path detection, and ToS handling that only runs for
    standard conda. Package-managed micromamba environments now use a
    no-space cache root to avoid `micromamba run -p` path parsing
    failures.
    [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
    now also stops early if reticulate has already initialized a
    different Python executable, avoiding unsafe in-session switches
    between conda and micromamba. The Python stack now also pins
    `setuptools < 81` so legacy packages such as `trimap` can still
    import `pkg_resources`.
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    The default `modules = NULL` environment no longer includes `scomm`;
    install it explicitly with `modules = "scomm"` because the
    TensorFlow/scOMM stack conflicts with the default JAX/scVI stack
    through incompatible `ml-dtypes` requirements.
  - [`CheckDataList()`](https://mengxu98.github.io/scop/reference/CheckDataList.md)
    /
    [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md):
    Removed an unnecessary all-feature
    [`ScaleData()`](https://mengxu98.github.io/scop/reference/ScaleData.md)
    call immediately after `LogNormalize`. Downstream workflows still
    scale the selected HVFs before PCA, but large raw-count inputs no
    longer create a dense all-gene `scale.data` intermediate during
    checking.
  - [`srt_reorder()`](https://mengxu98.github.io/scop/reference/srt_reorder.md):
    Replaced per-cluster repeated sparse subsetting and
    [`rowMeans()`](https://rdrr.io/pkg/Matrix/man/colSums-methods.html)
    with a single sparse group-membership matrix multiplication,
    speeding the cluster-reordering step used by
    [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md)
    while preserving average expression values.
  - [`FindExpressedMarkers()`](https://mengxu98.github.io/scop/reference/FindExpressedMarkers.md):
    Added a sparse fold-change and detection-rate prefilter path for
    common sparse RNA inputs, avoiding dense all-feature expression
    materialization before marker filtering while preserving default
    results. Supplied features are now intersected with the active layer
    before dense fallback scoring, which also keeps `scale.data` and
    partial-feature layers from passing absent rows into fold-change
    calculation.
  - C++-accelerated backends are now the default where available:
    [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md),
    [`RunGSVA()`](https://mengxu98.github.io/scop/reference/RunGSVA.md),
    [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md),
    [`RunDynamicEnrichment()`](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md),
    and
    [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md)
    now prefer `backend = "cpp"` while retaining `backend = "r"` for
    exact legacy/package behavior.
    [`RunPermutation()`](https://mengxu98.github.io/scop/reference/RunPermutation.md)
    now uses its validated C++ implementation directly.
  - [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md)
    and
    [`RunGSVA()`](https://mengxu98.github.io/scop/reference/RunGSVA.md):
    Added a reusable C++ gene-set scoring backend for `ssGSEA`,
    `zscore`, and `plage`;
    [`RunGSVA()`](https://mengxu98.github.io/scop/reference/RunGSVA.md)
    can now use `backend = "cpp"` for `method = "ssgsea"`,
    `method = "zscore"`, `method = "plage"`, and Gaussian- or
    Poisson-kernel `method = "gsva"`. PLAGE scores are oriented by the
    gene set mean z-score to avoid backend-dependent sign flips.
  - [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md)
    and
    [`RunGSVA()`](https://mengxu98.github.io/scop/reference/RunGSVA.md):
    Added `cpp_chunk_size` for the C++ GSVA kernel paths to reduce peak
    dense intermediate memory on large cell counts; `NULL` now
    auto-selects a chunk size for large matrices.
  - [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md):
    Added an experimental `backend = "cpp"` ORA path using a C++
    hypergeometric implementation for faster enrichment tables while
    keeping `backend = "r"` available as the clusterProfiler-compatible
    path.
  - [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md):
    Added experimental `backend = "cpp"` support for Seurat module
    scoring by keeping control-gene sampling in R and moving sparse mean
    calculations to C++ code.
  - [`RunPermutation()`](https://mengxu98.github.io/scop/reference/RunPermutation.md):
    Uses the C++ permutation and bootstrap loops directly after
    validating they match the legacy R calculation for observed
    fractions and log2 fold-differences while running substantially
    faster.
  - [`RunUMAP2()`](https://mengxu98.github.io/scop/reference/RunUMAP2.md):
    Added an internal C++ sparse column top-k helper for Graph inputs to
    speed extraction of precomputed neighbor indices/connectivities
    before calling `uwot`.
  - [`RunKNNMap()`](https://mengxu98.github.io/scop/reference/RunKNNMap.md)
    and
    [`RunKNNPredict()`](https://mengxu98.github.io/scop/reference/RunKNNPredict.md):
    Added an internal C++ dense column top-k helper for the raw KNN
    fallback to speed nearest-neighbor extraction from precomputed
    distance matrices.
  - [`PseudotimeProjectionPlot()`](https://mengxu98.github.io/scop/reference/PseudotimeProjectionPlot.md):
    Reused the internal dense top-k helper for pseudotime KNN/gradient
    neighbor extraction from distance matrices.
  - Mapping/integration metrics: Added an internal C++ contingency-table
    helper for `metric_accuracy()`, `metric_macro_f1()`,
    `metric_purity()`, `metric_nmi()`, `metric_ari()`,
    `metric_weighted_recall()`, and `collect_mapping_metrics()`.
  - [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md):
    Added a C++ backend for `method = "GSVA"` that keeps the
    Poisson-kernel GSVA scoring path but accelerates repeated
    count-kernel calculations.
  - [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md)
    and
    [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md):
    Added an experimental C++ backend for AUCell scoring selected only
    through `backend = "cpp"`, while keeping the original R/AUCell
    implementation available via `backend = "r"` for exact package
    output.
    [`RunDynamicEnrichment()`](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)
    now passes this backend through to AUCell-based scoring.
  - `CellScoring(method = "AUCell")`: Fixed score-name assignment when
    only one feature list is scored, avoiding vector dropping before
    metadata column names are applied.
  - `RunMetabolism(method = "VISION")`: Fixed signature construction so
    in-memory gene sets are passed as `VISION::Signature` objects
    instead of being interpreted as file paths.
  - [`DynamicHeatmap()`](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md)
    / custom enrichment workflows: Fixed incorrect term label display
    when user-provided `TERM2GENE` / `TERM2NAME` use non-standard column
    names such as `term` / `name`. Custom term annotations are now
    normalized consistently so that heatmap term labels show readable
    term names instead of fallback IDs (e.g. GO IDs). Related issue
    [\#160](https://github.com/mengxu98/scop/issues/160)
    ([@Pineapple-wen6](https://github.com/Pineapple-wen6),
    [@mengxu98](https://github.com/mengxu98)).
  - Refactored repeated custom database input handling by introducing
    shared internal helpers for normalizing and assembling user-provided
    `TERM2GENE` / `TERM2NAME`, and applied them across
    [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
    [`RunDynamicEnrichment()`](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md),
    [`RunGSEA()`](https://mengxu98.github.io/scop/reference/RunGSEA.md),
    [`RunGSVA()`](https://mengxu98.github.io/scop/reference/RunGSVA.md),
    and
    [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md)
    to keep behavior consistent.
  - [`RunIntegration()`](https://mengxu98.github.io/scop/reference/RunIntegration.md):
    Improved `ChromatinAssay` handling by standardizing ATAC reductions
    after integration, avoiding non-interactive small-batch blocking,
    clipping `TFIDF/rlsi` anchor dims to available cells, auto-switching
    `Harmony5` to legacy `Harmony`, and rejecting unsupported `Seurat` /
    `RPCA` ATAC paths with explicit messages.

## scop 0.8.7

- **feat**:
  - Added
    [`RunDimsEstimate()`](https://mengxu98.github.io/scop/reference/RunDimsEstimate.md)
    and
    [`DimsEstimatePlot()`](https://mengxu98.github.io/scop/reference/DimsEstimatePlot.md)
    for intrinsic dimensionality estimation from reductions, integrated
    into
    [`RunDimsReduction()`](https://mengxu98.github.io/scop/reference/RunDimsReduction.md)
    and
    [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md)
    to automatically select useful dimensions when
    `linear_reduction_dims_use = NULL`.
  - Renamed `RunDimReduction()` to
    [`RunDimsReduction()`](https://mengxu98.github.io/scop/reference/RunDimsReduction.md)
    and updated downstream callers/documentation accordingly.
  - [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md):
    When `linear_reduction_dims_use = NULL`, now uses estimated
    dimensions stored in the reduction (via
    [`RunDimsEstimate()`](https://mengxu98.github.io/scop/reference/RunDimsEstimate.md))
    when available, with fallback to the first 50 dimensions.
  - Added Seurat v5 integration methods:
    [`CCA_integrate()`](https://mengxu98.github.io/scop/reference/CCA_integrate.md),
    [`RPCA_integrate()`](https://mengxu98.github.io/scop/reference/RPCA_integrate.md),
    [`fastMNN5_integrate()`](https://mengxu98.github.io/scop/reference/fastMNN5_integrate.md),
    [`Harmony5_integrate()`](https://mengxu98.github.io/scop/reference/Harmony5_integrate.md),
    and
    [`scVI5_integrate()`](https://mengxu98.github.io/scop/reference/scVI5_integrate.md)
    via
    [`Seurat::IntegrateLayers()`](https://satijalab.org/seurat/reference/IntegrateLayers.html),
    and exposed them through
    [`RunIntegration()`](https://mengxu98.github.io/scop/reference/RunIntegration.md).
  - Added
    [`Coralysis_integrate()`](https://mengxu98.github.io/scop/reference/Coralysis_integrate.md)
    and exposed Coralysis through
    [`RunIntegration()`](https://mengxu98.github.io/scop/reference/RunIntegration.md),
    following the official Seurat v5 compatible workflow via
    `SingleCellExperiment`.
  - Added
    [`h5ad_to_srt()`](https://mengxu98.github.io/scop/reference/h5ad_to_srt.md)
    for reading `.h5ad` files directly into `Seurat` objects via
    `scanpy.read_h5ad()`, with automatic CSR/float64 coercion to avoid
    reticulate conversion issues. Layers that fail conversion are
    gracefully skipped and reported.
  - [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md):
    Improved robustness of layer conversion — each layer is now wrapped
    in [`tryCatch()`](https://rdrr.io/r/base/conditions.html) so that
    individual failures no longer abort the entire conversion; skipped
    layers are reported as warnings.
  - [`RunCellChat()`](https://mengxu98.github.io/scop/reference/RunCellChat.md):
    Enhanced to support condition-specific analyses and pairwise merged
    comparisons via `group_column` and `group_cmp` parameters.
  - Added
    [`RunCellphoneDB()`](https://mengxu98.github.io/scop/reference/RunCellphoneDB.md)
    for running CellphoneDB cell-cell communication analysis on a
    `Seurat` object through the official Python package, with support
    for species conversion and results stored in
    `srt@tools[["CellphoneDB"]]`.
  - Added
    [`RunNichenetr()`](https://mengxu98.github.io/scop/reference/RunNichenetr.md)
    and
    [`RunMultiNichenetr()`](https://mengxu98.github.io/scop/reference/RunMultiNichenetr.md)
    for running NicheNet and MultiNicheNet analysis on `Seurat` objects
    with standardized result storage.
  - Added unified cell-cell communication plotting functions
    [`CCCStatPlot()`](https://mengxu98.github.io/scop/reference/CCCStatPlot.md),
    [`CCCHeatmap()`](https://mengxu98.github.io/scop/reference/CCCHeatmap.md),
    and
    [`CCCNetworkPlot()`](https://mengxu98.github.io/scop/reference/CCCNetworkPlot.md)
    for CellChat, CellphoneDB, NicheNet, and MultiNicheNet results.
  - Remove `CellChatPlot()`.
- **fix**:
  - [`RunUMAP2()`](https://mengxu98.github.io/scop/reference/RunUMAP2.md):
    Fixed reduction lookup to check existing reduction names before
    falling back to
    [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md),
    avoiding errors when the exact reduction name is already present.

## scop 0.8.6

- **feat**:
  - Added
    [`RunGSVA()`](https://mengxu98.github.io/scop/reference/RunGSVA.md)
    for gene set variation analysis and
    [`GSVAPlot()`](https://mengxu98.github.io/scop/reference/GSVAPlot.md)
    for visualization of [GSVA](https://github.com/rcastelo/GSVA)
    results. Related issue
    [\#146](https://github.com/mengxu98/scop/issues/146)
    ([@mengxu98](https://github.com/mengxu98)).
  - [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md)
    - Added
      [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md)
      for single-cell metabolism pathway scoring with support for
      `AUCell`, `GSVA`, `ssGSEA`, and optional `VISION`, and
      [`MetabolismPlot()`](https://mengxu98.github.io/scop/reference/MetabolismPlot.md)
      for visualization
    - [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md)
      now uses [scMetabolism](https://github.com/wu-yc/scMetabolism)
      KEGG / Reactome metabolism pathway definitions to identify
      metabolism-related terms, then rebuilds pathway gene sets from
      updated
      [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md)
      annotations, enabling cached database updates and cross-species
      metabolism scoring. Related issue
      [\#146](https://github.com/mengxu98/scop/issues/146)
      ([@mengxu98](https://github.com/mengxu98)).
  - Added
    [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md)
    and integrated optional `decontX` ambient RNA decontamination into
    [`RunCellQC()`](https://mengxu98.github.io/scop/reference/RunCellQC.md),
    including contamination metadata, optional decontaminated assay
    output, and threshold-based QC filtering. Related issue
    [\#147](https://github.com/mengxu98/scop/issues/147)
    ([@mengxu98](https://github.com/mengxu98)).
- **fix**:
  - [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md):
    Fixed duplicated X/Y axis titles and main title when `stack = TRUE`
    with `theme_args` or `title`. Stack assembly now strips axis and
    plot titles from non-edge panels and draws a single shared
    ylab/title via gtable; theme styling from `theme_args`
    (e.g. `axis.title.y`, `title`/`plot.title`) is applied to these
    shared grobs. Related issue
    [\#145](https://github.com/mengxu98/scop/issues/145)
    ([@PanSX-Dr](https://github.com/PanSX-Dr)).
- **docs**:
  - Updated the README badge/logo.
- **data**:
  - Removed the example datasets `ifnb_sub`, `ref_scHCL`, and
    `ref_scZCL`.

## scop 0.8.5

- **feat**:
  - Default palette changed from `"Paired"` to `"Chinese"`. See
    [`thisplot::ChineseColors()`](https://mengxu98.github.io/thisplot/reference/ChineseColors.html)
    for details.
  - Unified the `cores` parameter across multiple functions
    ([`DynamicHeatmap()`](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md),
    [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md),
    [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md),
    and `heatmap_enrichment()`), ensuring `cores` is correctly threaded
    through to
    [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md).
  - [`RunMonocle2()`](https://mengxu98.github.io/scop/reference/RunMonocle2.md)
    and
    [`RunMonocle3()`](https://mengxu98.github.io/scop/reference/RunMonocle3.md):
    Added `xlab` and `ylab` parameters, passed through to internal
    [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
    and
    [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
    calls.
  - [`ListDB()`](https://mengxu98.github.io/scop/reference/ListDB.md):
    Now supports multiple species in the `species` parameter
    simultaneously and adds `Species` and `DB` columns to the output
    data frame for clearer identification.
  - Moved `is_outlier` to
    [`thisutils::is_outlier()`](https://mengxu98.github.io/thisutils/reference/is_outlier.html).
  - `configure_apple_silicon_env()` (*Python* function): Added OpenMP
    compatibility handling on macOS arm64 by prepending environment
    `lib` paths to `DYLD_FALLBACK_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` and
    preloading `libomp.dylib`/`libiomp5.dylib` via
    `ctypes.CDLL(..., RTLD_GLOBAL)` to reduce `scanpy`/`python-igraph`
    import conflicts.
  - [`scVI_integrate()`](https://mengxu98.github.io/scop/reference/scVI_integrate.md):
    Unified parameter naming by renaming `num_threads` to `cores` for
    consistency across integration functions.
  - Added `find_neighbors_and_clusters()` and
    `run_nonlinear_reduction()` helper functions to unify operations
    across integration functions and reduce redundant code.
- **fix**:
  - [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md):
    - TF database now uses
      [AnimalTFDB4](https://github.com/mengxu98/datasets/tree/main/AnimalTFDB4)
      as the data source.
    - MP database - the Web Archive URL year is now dynamically
      determined from the current date, and the archive snapshot date
      for all MP-related file downloads (`VOC_MammalianPhenotype.rpt`,
      `MGI_Gene_Model_Coord.rpt`, `MGI_GenePheno.rpt`) is extracted from
      the server index page to ensure consistent versioning.
    - hTFtarget data download URL has been changed from
      `"http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"`
      to
      `"https://guolab.wchscu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"`.
    - CSPA data download URL has been changed from
      `"https://wlab.ethz.ch/cspa/data/S1_File.xlsx"` to
      `"https://raw.githubusercontent.com/mengxu98/CSPA/main/S1_File.xlsx"`.

    Related issus [\#76](https://github.com/mengxu98/scop/issues/76)
    ([@hwa2Hu](https://github.com/hwa2Hu)),
    [\#139](https://github.com/mengxu98/scop/issues/139)
    ([@pengding774-dot](https://github.com/pengding774-dot)),
    [\#140](https://github.com/mengxu98/scop/issues/140)
    ([@mengxu98](https://github.com/mengxu98)).
  - [`LIGER_integrate()`](https://mengxu98.github.io/scop/reference/LIGER_integrate.md)
    - Migrated to the `rliger` 2.x workflow
      ([`rliger::runIntegration()`](https://welch-lab.github.io/liger/reference/runIntegration.html) +
      [`rliger::quantileNorm()`](https://welch-lab.github.io/liger/reference/quantileNorm.html)
      on `Seurat` object) and now prepares/uses the `ligerScaleData`
      layer via
      [`rliger::scaleNotCenter()`](https://welch-lab.github.io/liger/reference/scaleNotCenter.html)
      before integration.
    - Updated argument naming from `LIGER_dims_use` to `liger_dims_use`,
      and removed legacy quantile-normalization parameter compatibility
      mapping (`ref_dataset`), keeping `reference` as the supported
      interface.
  - Optimized the installation of some *Python* packages on Apple
    Silicon devices.
- **docs**:
  - Updated README example code and visualizations. After regenerating
    figures, the package size was reduced by ~3 MB.

## scop 0.8.4

- **fix**:
  - [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)
    / `ExpressionStatPlot()`: Fixed box and violin x-axis misalignment
    when `add_box = TRUE` with `split.by`. Groups with fewer than 2
    observations are now filtered before violin density estimation (with
    a warning), and the violin layer uses a consistent
    `position_dodge(width = 0.9)` to match the boxplot. Related issue
    [\#123](https://github.com/mengxu98/scop/issues/123)
    ([@oranges7](https://github.com/oranges7)).

## scop 0.8.3

- **feat**:
  - [`RunDynamicFeatures()`](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md):
    Added [PreTSA](https://github.com/haotian-zhuang/PreTSA/) method for
    dynamic feature fitting. The
    [PreTSA](https://github.com/haotian-zhuang/PreTSA/) algorithm in
    `scop` has been re implemented to support parallelization for higher
    performance. Original research: [PreTSA: computationally efficient
    modeling of temporal and spatial gene expression
    patterns](https://doi.org/10.1186/s13059-026-03994-3). Use
    `fit_method = "pretsa"` for B-spline-based piecewise truncated
    spline analysis; `fit_method = "gam"` (default) keeps generalized
    additive models. PreTSA supports `knot` (`0` or `"auto"`) and
    `max_knot_allowed` when `knot = "auto"`. Relate issue
    [\#133](https://github.com/mengxu98/scop/issues/133).
  - [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
    and
    [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md):
    Added `legend.title` parameter (default `NULL`) to control the
    legend title. When `NULL`, default titles are used (e.g. group name
    for `CellDimPlot`, feature or empty for `FeatureDimPlot`).
  - `ExpressionStatPlot()` and
    [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md):
    Added `legend.title` parameter (default `NULL`) for single-legend
    plots. When `NULL`, the default title (e.g. `keynm` or feature/group
    name) is used.
    [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)
    forwards `legend.title` to `ExpressionStatPlot()`.
  - Moved `StatPlot` function to
    [`thisplot::StatPlot`](https://mengxu98.github.io/thisplot/reference/StatPlot.html).
- **fix**:
  - Differential expression plots
    ([`VolcanoPlot()`](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
    [`DEtestManhattanPlot()`](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md),
    [`DEtestRingPlot()`](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)):
    added `only.pos = TRUE` for positive-only visualization.
    [`DEtestManhattanPlot()`](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md)
    now keeps the cell-type track centered at y = 0 and sizes the track
    from the nearest point distance around zero.
  - [`DynamicHeatmap()`](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md)
    / `heatmap_enrichment()`: Fixed incorrect `db` handling when using
    custom `TERM2GENE`/`TERM2NAME`. Enrichment results with
    `Database = "custom"` could be incorrectly filtered by default `db`
    values (e.g. `"GO_BP"`), causing false “No term enriched using the
    threshold” warnings even when enrichment succeeded. Relate issue
    [\#133](https://github.com/mengxu98/scop/issues/133)
    ([@1228849000](https://github.com/1228849000)).

## scop 0.8.2

- **feat**:
  - Add the
    [`DEtestPlot()`](https://mengxu98.github.io/scop/reference/DEtestPlot.md)
    function, which calls the original
    [`VolcanoPlot()`](https://mengxu98.github.io/scop/reference/VolcanoPlot.md)
    and adds two plot types, Manhattan and Ring, controlled by the
    `plot_type` parameter (`c("volcano", "manhattan", "ring")`, default
    `"volcano"`). Add standalone functions
    [`DEtestManhattanPlot()`](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md)
    and
    [`DEtestRingPlot()`](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)
    for direct use. Relate issue
    [\#121](https://github.com/mengxu98/scop/issues/121)
    ([@ericavalentini](https://github.com/ericavalentini)).
  - Differential expression visualization
    ([`DEtestPlot()`](https://mengxu98.github.io/scop/reference/DEtestPlot.md),
    [`VolcanoPlot()`](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
    [`DEtestManhattanPlot()`](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md),
    [`DEtestRingPlot()`](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)):
    added `res` parameter to accept existing DE results (data.frame).
    When `res` is provided, `srt` is ignored. Data processing
    supports: (1) `group1` or `cluster` column for grouped plots; (2) no
    grouping column for a single panel; (3) gene names from row names
    when `gene` column is missing. Relate issue
    [\#129](https://github.com/mengxu98/scop/issues/129)
    ([@mengxu98](https://github.com/mengxu98)).
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    Update `version` parameter to specify the *Python* version of the
    conda environment. Default is `"3.11-1"` on *Windows* and `"3.10-1"`
    on *macOS* and *Unix*. Relate issue
    [\#103](https://github.com/mengxu98/scop/issues/103)
    ([@PanSX-Dr](https://github.com/PanSX-Dr)).
  - [`RunNMF()`](https://mengxu98.github.io/scop/reference/RunNMF.md):
    Add the `cores` parameter for
    [`RunNMF()`](https://mengxu98.github.io/scop/reference/RunNMF.md)
    and optimize the printed message.
  - Removed the setting in *Python* functions that prevents drawing
    functions from causing *R* crashes.
- **fix**:
  - [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md):
    Fixed unused `plot_format` parameter error. The parameter is now
    properly excluded from arguments passed to Python functions. Relate
    issue [\#126](https://github.com/mengxu98/scop/issues/126)
    ([@christinejay990202-dev](https://github.com/christinejay990202-dev)).

## scop 0.8.1

- **fix**:
  - [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md):
    Fixed `group_palcolor` when a named vector is passed: the function
    previously used only the first color for all groups because
    `group_palcolor[[1]]` on a vector returns a single element. Now when
    `group.by` has length 1, a vector is automatically wrapped as a
    list; when `within_groups = TRUE`, `group_palcolor` is expanded in
    line with `group_palette`.
  - [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md):
    Same `group_palcolor` fix as
    [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md):
    support for named-vector input and correct expansion when
    `within_groups = TRUE`.
- **docs**:
  - `group_by` modified to `group.by` to unify documentation, involving
    functions:
    [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md),
    [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md),
    [`VolcanoPlot()`](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
    [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
    [`EnrichmentPlot()`](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md),
    [`RunGSEA()`](https://mengxu98.github.io/scop/reference/RunGSEA.md),
    [`GSEAPlot()`](https://mengxu98.github.io/scop/reference/GSEAPlot.md),
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md),
    [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md),
    [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md),
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
    releated to [\#120](https://github.com/mengxu98/scop/issues/120).

## scop 0.8.0

- **feat**:
  - [`RunMonocle2()`](https://mengxu98.github.io/scop/reference/RunMonocle2.md):
    New function for performing
    [Monocle2](https://github.com/mengxu98/monocle) trajectory analysis
    with support for various dimensionality reduction methods (DDRTree,
    ICA, tSNE, SimplePPT, L1-graph, SGL-tree). Uses the fixed version of
    monocle2 from
    [mengxu98/monocle](https://github.com/mengxu98/monocle).
  - [`RunMonocle3()`](https://mengxu98.github.io/scop/reference/RunMonocle3.md):
    New function for performing
    [Monocle3](https://github.com/cole-trapnell-lab/monocle3) trajectory
    analysis with support for cell ordering, trajectory learning, and
    pseudotime computation.
  - [`RunCytoTRACE()`](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md):
    New `scop` implementation for running CytoTRACE 2 analysis to
    predict cellular potency scores and categories (Differentiated,
    Unipotent, Oligopotent, Multipotent, Pluripotent, Totipotent) with
    support for human and mouse species.
  - [`CytoTRACEPlot()`](https://mengxu98.github.io/scop/reference/CytoTRACEPlot.md):
    New function for visualizing CytoTRACE 2 analysis results.

## scop 0.7.9

- **docs**:
  - Optimized `@inheritParams` usage to reduce redundant parameter
    definitions.
  - Unified parameter naming: renamed `group_by` to `group.by` for
    consistency with Seurat conventions in
    [`RunMonocle3()`](https://mengxu98.github.io/scop/reference/RunMonocle3.md),
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md),
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
    [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md),
    [`CellCorHeatmap()`](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md),
    and
    [`RunDynamicFeatures()`](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md).

## scop 0.7.8

- **fix**:
  - [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md):
    Fixed unused `plot_format` parameter error. The parameter is now
    properly excluded from arguments passed to Python functions. Relate
    issue [\#114](https://github.com/mengxu98/scop/issues/114)
    ([@Moonerss](https://github.com/Moonerss)).
  - [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md):
    Fixed PAGA computation error by replacing `scv.tl.paga` with
    `sc.tl.paga` (scanpy implementation) for better stability. The
    function now uses the same PAGA implementation as
    [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md)
    function.
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md):
    Fixed GPCCA Schur decomposition error by adding fallback mechanism.
    When `brandts` method fails with “subspace_angles” error, the
    function automatically tries `krylov` method. If both methods fail,
    it automatically switches to CFLARE estimator for more robust
    computation.
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md):
    Fixed `recover_dynamics` error by ensuring `velocity_graph` and
    `velocity_graph_neg` are properly set before calling
    `scv.tl.recover_dynamics()` for latent time computation.
- **feat**:
  - [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md):
    Removed automatic removal of “X\_” prefix from dimensionality
    reduction names in `obsm` keys. The function now preserves original
    reduction names as they are stored in AnnData objects.
- **data**:
  - Reducing the size of `pancreas_sub` example dataset.

## scop 0.7.7

- **feat**:
  - [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md):
    Enhanced to support multiple AnnData object types including *Python*
    AnnData objects (from scanpy/reticulate), R6 AnnData objects from
    the `anndata` package (AnnDataR6), and R6 AnnData objects from the
    `anndataR` package (InMemoryAnnData). Added internal helper
    functions `get_adata_element()` and `get_adata_names()` for better
    compatibility. Relate issues
    [\#67](https://github.com/mengxu98/scop/issues/67)
    ([@lisch7](https://github.com/lisch7)),
    [\#91](https://github.com/mengxu98/scop/issues/91)
    ([@mengxu98](https://github.com/mengxu98)) and
    [commit91#issuecomment](https://github.com/mengxu98/scop/issues/91#issuecomment-3659404993).
- **fix**:
  - [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md):
    Fixed error when comparing one cluster against multiple clusters
    using `group1` and `group2` parameters. Relate issue
    [\#111](https://github.com/mengxu98/scop/issues/111)
    ([@zhaoxiaoyan9225](https://github.com/zhaoxiaoyan9225)).
  - [`AnnotateFeatures()`](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md):
    Fixed bug where the function would fail when processing GTF file
    annotations due to column name matching issues during data naming.
    The function now correctly handles column name intersections when
    merging annotation data.

## scop 0.7.6

- **feat**:
  - [`RunDM()`](https://mengxu98.github.io/scop/reference/RunDM.md):
    Added automatic PCA-based dimensionality reduction when using many
    features (\>1000) to speed up diffusion map computation. The `npcs`
    parameter can be used to control the number of principal components
    used for pre-processing.

## scop 0.7.5

- **fix**:
  - [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md):
    Fixed bug where the function failed to build results. Relate issue
    [\#98](https://github.com/mengxu98/scop/issues/98)
    ([@SuperrNaruto](https://github.com/SuperrNaruto)).
- **feat**:
  - [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md):
    Fixed compatibility issue with
    [SeuratObject](https://satijalab.github.io/seurat-object/) 5.0.0+ by
    replacing deprecated `Assays()` `slot` argument with `LayerData()`.
    Relate issue [\#100](https://github.com/mengxu98/scop/issues/100)
    ([@mattizecos](https://github.com/mattizecos)).
  - [`RunDM()`](https://mengxu98.github.io/scop/reference/RunDM.md):
    Added automatic PCA-based dimensionality reduction when using many
    features (\>1000) to speed up diffusion map computation. The `npcs`
    parameter can be used to control the number of principal components
    used for pre-processing.

## scop 0.7.3

- **feat**:
  - [`RunCellTypist()`](https://mengxu98.github.io/scop/reference/RunCellTypist.md):
    New function for cell type annotation using the CellTypist method.
  - [`CellTypistModels()`](https://mengxu98.github.io/scop/reference/CellTypistModels.md):
    New function for downloading and managing CellTypist pre-trained
    models.

## scop 0.7.2

- **feat**:
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md):
    Performance optimizations and code improvements.

## scop 0.7.1

- **fix**:
  - [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md):
    Fixed issue where NA values appeared in labels. Relate issue
    [\#93](https://github.com/mengxu98/scop/issues/93)
    ([@12345nkjil](https://github.com/12345nkjil)).

## scop 0.7.0

- **feat**:
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    Integrated `uv` as the primary *Python* package installer for
    improved installation speed.
  - [`check_python()`](https://mengxu98.github.io/scop/reference/check_python.md):
    Now uses `uv` as the primary installation tool with `pip` as
    fallback, significantly improving package installation speed.
  - Added `find_uv()` and `install_uv()` internal functions for managing
    `uv` package manager installation and detection.

## scop 0.6.6

- **docs**:
  - Unified documentation format across all R functions:
    - Standardized return value tags: Changed all `@returns` to
      `@return` for consistency.
    - Unified parameter documentation: Replaced all `\code{value}` with
      Markdown backticks `` `value` `` format.
    - Standardized default value descriptions.
    - Added `@md` tags: Added `@md` tags to all functions using Markdown
      syntax in documentation.
    - Enhanced cross-references: Added `@seealso` links to related
      functions where appropriate.

## scop 0.6.5

- **feat**:
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    - Added comprehensive environment variable configuration to prevent
      crashes when calling *Python* functions, including setting thread
      limits for OMP, OPENBLAS, MKL, NUMBA, and other libraries. This
      improves stability on all platforms, especially Apple silicon
      Macs.
    - Added `accept_conda_tos()` function to automatically accept conda
      Terms of Service for required channels, improving the conda
      environment setup process.
    - Fixed conda Terms of Service acceptance issue in
      [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md).
      The function now automatically accepts conda Terms of Service for
      required channels, eliminating the need for manual acceptance.
      This addresses the issue reported in
      [\#85](https://github.com/mengxu98/scop/issues/85).
  - Multiple Python-based functions (`RunPAGA`, `RunSCVELO`,
    `RunPalantir`, `RunCellRank`, `RunWOT`, `RunPHATE`, `RunPaCMAP`,
    `RunTriMap`): Enhanced message formatting and code improvements.
  - [`PrepareSCExplorer()`](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md):
    Fixed package version dependency issues with `shiny` and `bslib`
    compatibility. The function now properly handles `bslib` theme
    configuration to work with both `shiny` 1.6.0 and 1.7.0+, addressing
    compatibility errors reported in
    [\#87](https://github.com/mengxu98/scop/issues/87).
- **fix**:
  - Improved code formatting and consistency across multiple functions.
  - Enhanced Python functions in `inst/python/functions.py` with better
    error handling and message formatting.
- **docs**:
  - Updated documentation for multiple functions to reflect code
    improvements.

## scop 0.6.2

- **feat**:
  - `CellChatPlot()`: Adjusted the size of saved figures for better file
    size optimization.
- **docs**:
  - Updated README.md to remove references to Monocle2 and Monocle3
    (deprecated functions).

## scop 0.6.1

- **feat**:
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    Improved message formatting and simplified log output for better
    user experience.
  - Added `get_conda_envs_dir()` helper function to centralize conda
    environment directory retrieval.
  - [`RunIntegration()`](https://mengxu98.github.io/scop/reference/RunIntegration.md):
    Enhanced `integration_method` parameter definition with explicit
    method list for better code clarity.
- **fix**:
  - Moved `exist_python_pkgs()` function to `check_package.R` for better
    code organization.
  - Replaced direct `conda_info()$envs_dirs[1]` calls with
    `get_conda_envs_dir()` helper function for consistency.
  - [`RunSCExplorer()`](https://mengxu98.github.io/scop/reference/RunSCExplorer.md):
    Updated to use
    [`thisplot::palette_list`](https://mengxu98.github.io/thisplot/reference/palette_list.html)
    and
    [`thisplot::slim_data()`](https://mengxu98.github.io/thisplot/reference/slim_data.html)
    instead of `scop::palette_list` and `scop::slim_data()`.
  - Added `thisplot` to dependency checks in
    [`RunSCExplorer()`](https://mengxu98.github.io/scop/reference/RunSCExplorer.md).
- **docs**:
  - Updated documentation across multiple functions.

## scop 0.6.0

- **feat**:
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    Enhanced with environment caching mechanism to avoid redundant
    environment preparation. Improved message formatting and error
    handling.
  - Python-based functions
    ([`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md),
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md),
    [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md),
    [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md))
    now automatically call
    [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
    internally, eliminating the need for users to manually prepare the
    Python environment before using these functions.
  - `cluster_within_group2()`: New function for clustering within
    groups.
  - Multiple plotting functions: Replaced `geom_sankey()` with
    `ggsankey::geom_sankey()` for better Sankey diagram support.
  - Multiple functions: Replaced `:::` operator with
    `get_namespace_fun()` for safer namespace access.
- **fix**:
  - Removed `RunMonocle()` function and related documentation
    (`RunMonocle2.Rd`, `RunMonocle3.Rd`).
  - Removed `projection_functions.R` file (functions moved to other
    locations).
  - Replaced custom theme functions with
    [`thisplot::theme_this()`](https://mengxu98.github.io/thisplot/reference/theme_this.html)
    (exported as
    [`theme_scop()`](https://mengxu98.github.io/scop/reference/theme_scop.md)).
  - Replaced direct `log_message()` calls with
    [`thisutils::log_message()`](https://mengxu98.github.io/thisutils/reference/log_message.html)
    for consistency.
  - Removed `palette_list` data object.
- **deps**:
  - Moved `cli` from `Suggests` to `Imports` for better message
    formatting support.
  - Added `thisplot` to `Imports` for theme and utility functions.
  - Added `ggsankey` to `Suggests` for Sankey diagram support.
  - Added remote dependencies: `theislab/destiny`, `mengxu98/thisplot`.
  - Removed unused dependencies: `Biobase`, `BiocGenerics`,
    `concaveman`, `DDRTree`, `glmGamPoi`, `hexbin`, `monocle`, `png`,
    `ragg`, `tidyr`.
- **docs**:
  - Updated documentation across multiple functions to reflect code
    refactoring.
  - Improved code organization and maintainability.

## scop 0.5.5

- **fix**:
  - Fixed
    [`VelocityPlot()`](https://mengxu98.github.io/scop/reference/VelocityPlot.md)
    function error in `plot_type = "grid"` mode: replaced vectorized
    arrow length with fixed-length arrows (using mean length) to resolve
    [`vapply()`](https://rdrr.io/r/base/lapply.html) error that occurred
    when [`grid::arrow()`](https://rdrr.io/r/grid/arrow.html) received a
    vector instead of a single value, see
    [\#74](https://github.com/mengxu98/scop/issues/74).

## scop 0.5.4

- **fix**:
  - Fixed parameter name error in
    [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
    function calls: changed `data` parameter to `object` in
    [`RunKNNMap()`](https://mengxu98.github.io/scop/reference/RunKNNMap.md),
    [`RunSymphonyMap()`](https://mengxu98.github.io/scop/reference/RunSymphonyMap.md),
    [`RunScmap()`](https://mengxu98.github.io/scop/reference/RunScmap.md),
    [`RunPCAMap()`](https://mengxu98.github.io/scop/reference/RunPCAMap.md),
    and
    [`RunSingleR()`](https://mengxu98.github.io/scop/reference/RunSingleR.md)
    functions, see [\#68](https://github.com/mengxu98/scop/issues/68).
  - Fixed `SingleCellExperiment` object creation in
    [`RunScmap()`](https://mengxu98.github.io/scop/reference/RunScmap.md)
    and
    [`RunSingleR()`](https://mengxu98.github.io/scop/reference/RunSingleR.md)
    functions: changed from coercing `SummarizedExperiment` to directly
    constructing `SingleCellExperiment` objects.

## scop 0.5.3

- **feat**:
  - [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md):
    Changed default `Ensembl_version` parameter from `103` to `NULL` for
    more flexible version handling.
  - Added *Python* version `log_message()` for Python-based functions
    ([`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
    [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md),
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md),
    [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md),
    [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md))
    and added `verbose` parameter inheritance and improved message
    formatting through cli formatting.
- **fix**:
  - Delete `harmonizomeapi.py` file.
  - Move `scop_analysis.py` into a single `functions.py` file in
    `inst/python/` for better code organization and maintainability.
- **docs**:
  - Improved parameter documentation consistency.

## scop 0.5.1

- **docs**:
  - Improved reference formatting and consistency across multiple
    functions.
  - Enhanced documentation clarity and readability.

## scop 0.5.0

- **feat**:
  - [`RunCellChat()`](https://mengxu98.github.io/scop/reference/RunCellChat.md):
    New function to perform CellChat analysis for investigating
    cell-to-cell communication with support for human, mouse, and
    zebrafish species.
  - `CellChatPlot()`: New function to visualize CellChat analysis
    results with various plot types and customization options.
  - Multiple integration functions: Improved error messages and message
    formatting for better user experience.
- **deps**:
  - Added [CellChat](https://github.com/jinworks/CellChat) package
    dependency with remote repository `jinworks/CellChat`.
- **docs**:
  - Updated README.md with improved code formatting and examples.
  - Enhanced documentation for cell communication analysis functions.
  - Improved error messages and user guidance across integration
    functions.
- **fix**:
  - Removed some example figures to optimize package installation size.

## scop 0.4.0

- **feat**:
  - [`RunProportionTest()`](https://mengxu98.github.io/scop/reference/RunProportionTest.md):
    New function to perform Monte-carlo permutation test for quantifying
    cell proportion differences between conditions.
  - [`ProportionTestPlot()`](https://mengxu98.github.io/scop/reference/ProportionTestPlot.md):
    New function to generate proportion test plots with customizable
    significance thresholds and visualization options.
  - Multiple *Python*-based functions: add `\dontrun{}` blocks for
    Github workfolw checking.
- **docs**:
  - Added comprehensive documentation for new proportion testing
    functions.
  - Enhanced example documentation across multiple functions.
  - Updated package documentation and examples.

## scop 0.3.4

- **docs**:
  - Updated workflow examples and function documentation.

## scop 0.3.3

- **feat**:
  - Multiple functions: Improved parameter documentation formatting and
    consistency across the package.

## scop 0.3.2

- **feat**:
  - [`GetFeaturesData()`](https://mengxu98.github.io/scop/reference/GetFeaturesData.md)
    and
    [`AddFeaturesData()`](https://mengxu98.github.io/scop/reference/AddFeaturesData.md):
    Enhanced argument clarity, added input validation, and standardized
    return values for `Seurat`, `Assay`, and `Assay5` objects.
  - [`CellCorHeatmap()`](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md):
    - Renamed parameters: `query_cell_annotation` → `query_annotation`,
      `ref_cell_annotation` → `ref_annotation`.
    - Improved error message formatting through cli formatting.
    - Simplified variable assignments and improved readability.
- **docs**:
  - Comprehensive documentation updates across multiple functions
    including `AnnotateFeatures`, `CellDimPlot`, `CellStatPlot`,
    `FeatureStatPlot`, `GroupHeatmap`, `RunCellQC`, and others.
  - Improved parameter descriptions and function clarity.

## scop 0.3.1

- **feat**:
  - [`EnrichmentPlot()`](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md)
    and
    [`GSEAPlot()`](https://mengxu98.github.io/scop/reference/GSEAPlot.md):
    Removed conditional font face styling (`face = ifelse()` logic) for
    better text rendering consistency. Set the default value of
    `lineheight` from `0.5` to `0.7`.
  - Updated `check_r()` function for improved package checking
    functionality.
  - Updated reexports functionality.
- **docs**:
  - Updated documentation formatting and consistency.

## scop 0.3.0

- **feat**:
  - Fixed `segmentation faults` and `R crashes` on *M-series* MacBook
    when running *Python* functions.
  - [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md):
    Enhanced with *M-series* MacBook detection and automatic environment
    configuration.
  - [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md):
    Added ARM64-specific optimizations to prevent crashes and ensure
    stable execution.
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md):
    Implemented *M-series* compatibility with proper NUMBA
    configuration.
  - [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md):
    Added ARM64 support with single-threaded execution mode.
  - [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md):
    Enhanced with *M-series* MacBook environment variable settings.
  - [`RunTriMap()`](https://mengxu98.github.io/scop/reference/RunTriMap.md):
    Added *M-series* MacBook compatibility for dimensionality reduction.
  - [`RunPaCMAP()`](https://mengxu98.github.io/scop/reference/RunPaCMAP.md):
    Implemented ARM64-specific environment configuration.
  - [`RunPHATE()`](https://mengxu98.github.io/scop/reference/RunPHATE.md):
    Added *M-series* MacBook support for non-linear dimensionality
    reduction.
  - [`RunCellQC()`](https://mengxu98.github.io/scop/reference/RunCellQC.md):
    Enhanced both scrublet and doubletdetection functions with ARM64
    compatibility.
- **docs**:
  - Updated function documentation to reflect *M-series* MacBook
    compatibility.
  - Added technical notes about ARM64 architecture considerations.

## scop 0.2.9

- **fix**:
  - Fix bug for
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md).
- **docs**:
  - Updated documentation for some functions.

## scop 0.2.7

- **feat**:
  - Added an internal function `.check_pkg_status()` to check if an *R*
    package is installed.
  - Update function
    [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
    to *S4* class function.
  - Update function
    [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md),
    make it more efficient.
- **data**:
  - Delete `lifemap` data, including: `lifemap_cell`,
    `lifemap_compartment` and `lifemap_organ`.
  - Reconstructed the sample data for both `panc8_sub` and
    `pancreas_sub`, retaining only the basic `Seurat` object.

## scop 0.2.6

- **feat**:
  - Added `remove_r()` function for easy remove *R* packages.
  - Rename function: `RemovePackages()` to
    [`remove_python()`](https://mengxu98.github.io/scop/reference/remove_python.md).
  - Removed other methods of installing *R* packages from the
    `check_r()` function, only retaining
    [pak::pak](https://pak.r-lib.org/reference/pak.html).
  - Delete useless import packages: `BBmisc`, `BiocManager`, `covr`,
    `devtools`, `promises` and `withr`.
  - Optimize the structure of `_pkgdown.yml` file.
- **docs**:
  - Updated documentation for some functions.

## scop 0.2.5

- **feat**:
  - Rename function: `palette_scop()` to `palette_colors()`.

## scop 0.2.4

- **feat**:
  - Rename functions: `check_srt_merge()` to
    [`CheckDataMerge()`](https://mengxu98.github.io/scop/reference/CheckDataMerge.md),
    `check_srt_list()` to `CheckDataList` and `check_data_type()` to
    `CheckDataType`.

## scop 0.2.2

- **feat**:
  - Replace all
    [`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html)
    with
    [`thisutils::parallelize_fun()`](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html).
- **fix**:
  - Fix bugs in
    [`RunSingleR()`](https://mengxu98.github.io/scop/reference/RunSingleR.md).

## scop 0.2.0

- **feat**:
  - Added
    [`remove_python()`](https://mengxu98.github.io/scop/reference/remove_python.md)
    function for easy remove *Python* packages.
- **fix**:
  - Corrected an issue in `py_to_r2()` function (intrinsic function),
    which ensures that Python-dependent functions like
    [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md)
    and
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md)
    function run correctly.

## scop 0.1.9

- **feat**:
  - Update
    [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md)
    and `AddModuleScore2()` functions. Now, new parameters `cores` and
    `verbose` have been added.
  - `AddModuleScore2()` function no longer uses the
    [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
    function to enable parallelization, but
    [`thisutils::parallelize_fun`](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html),
    and the `cores` parameter is used to control the number of cores in
    [`thisutils::parallelize_fun`](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html).

## scop 0.1.5

- **fix**:
  - Fix error for
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md)
    function, see [\#23](https://github.com/mengxu98/scop/issues/23).

## scop 0.1.4

- **feat**:
  - Update `.onAttach()`, now `.onAttach()` will print more information
    about *conda* and *Python*.
  - Update
    [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
    function for easy add or update a conda environments and install
    *Python* packages.
  - Added
    [`ListEnv()`](https://mengxu98.github.io/scop/reference/ListEnv.md)
    and
    [`RemoveEnv()`](https://mengxu98.github.io/scop/reference/RemoveEnv.md)
    functions for easy management of *conda* environment and *Python*
    packages.

## scop 0.1.3

- **feat**:
  - Added
    [`TACSPlot()`](https://mengxu98.github.io/scop/reference/TACSPlot.md)
    function for creating FACS-like plots. Please refer to [Kernfeld et
    al. paper](https://doi.org/10.1016/j.immuni.2018.04.015) and
    [Github](https://github.com/maehrlab/thymusatlastools2/blob/f8b51ad684d56b2eeda780787eb9ad4ff3003eef/R/data_handling_seurat.R#L271)
    for specific information.

## scop 0.0.9

- **fix**:
  - Fix a bug for
    [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
    function.
  - The default *Python* version is now set to `3.10-1`.

## scop 0.0.6

- **feat**:
  - Add
    [`GetAssayData5()`](https://mengxu98.github.io/scop/reference/GetAssayData5.md)
    function, a reimplementation of `GetAssayData()`, for compatibility
    with Seurat v5 `Assay` objects.
  - Updated
    [`GetFeaturesData()`](https://mengxu98.github.io/scop/reference/GetFeaturesData.md)
    and
    [`AddFeaturesData()`](https://mengxu98.github.io/scop/reference/AddFeaturesData.md)
    function to support retrieving and adding feature metadata.

## scop 0.0.5

- **data**:
  - Updated the `pancreas_sub` and `panc8_sub` test datasets to the
    Seurat v5 object format.

## scop 0.0.1

- **Initial version**
