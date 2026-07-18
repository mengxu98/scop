# Spatial method discovery and result inspection --------------------------------

.spatial_backend_cache <- new.env(parent = emptyenv())

spatial_registry_entry <- function(
  method,
  kind,
  task,
  implementation_files,
  tool_key = NA_character_,
  backend_id = "core",
  status = "stable",
  coordinate_space_current = "none",
  coordinate_space_target = coordinate_space_current,
  coordinate_requirement = "none",
  scalability = "standard",
  plot_function = NA_character_,
  documentation = method,
  backend_requirement = "all"
) {
  data.frame(
    method = method,
    kind = kind,
    task = task,
    implementation_files = implementation_files,
    tool_key = tool_key,
    backend_id = backend_id,
    status = status,
    coordinate_space_current = coordinate_space_current,
    coordinate_space_target = coordinate_space_target,
    coordinate_requirement = coordinate_requirement,
    scalability = scalability,
    plot_function = plot_function,
    documentation = documentation,
    backend_requirement = backend_requirement,
    stringsAsFactors = FALSE
  )
}

spatial_method_registry <- function() {
  entry <- spatial_registry_entry
  rows <- list(
    entry("SpatialCoordinates", "accessor", "data_access", "SpatialCore.R", backend_id = "core", coordinate_space_current = "mixed", coordinate_space_target = "mixed", coordinate_requirement = "backend_managed"),
    entry("GetSpatialResult", "accessor", "result_access", "SpatialRegistry.R", backend_id = "core", coordinate_space_current = "none"),
    entry("GetSpatialGraph", "accessor", "network", "RunSpatialNetwork.R", backend_id = "core", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("RunSpotQC", "analysis", "quality_control", "RunSpotQC.R", backend_id = "core", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "SpatialSpotPlot"),
    entry("RunSpaNorm", "analysis", "normalization", "RunSpaNorm.R", "SpaNorm", "spanorm", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),
    entry("RunSpatialQM", "analysis", "quality_control", "RunSpatialQM.R", "SpatialQM", "spatialqm", coordinate_space_current = "mixed", coordinate_space_target = "mixed", coordinate_requirement = "backend_managed", plot_function = "SpatialSpotPlot"),
    entry("RunSpotSweeper", "analysis", "quality_control", "RunSpotSweeper.R", "SpotSweeper", "spotsweeper", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),

    entry("srt_to_giotto", "bridge", "framework_bridge", "SpatialFrameworkConvert.R", backend_id = "giotto;giotto_class", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("giotto_to_srt", "bridge", "framework_bridge", "SpatialFrameworkConvert.R", backend_id = "giotto_class", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("srt_to_spata2", "bridge", "framework_bridge", "SpatialFrameworkConvert.R", backend_id = "spata2", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("spata2_to_srt", "bridge", "framework_bridge", "SpatialFrameworkConvert.R", backend_id = "spata2", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("srt_to_spe", "bridge", "framework_bridge", "SpatialDataConvert.R", backend_id = "core", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("spe_to_srt", "bridge", "framework_bridge", "SpatialDataConvert.R", backend_id = "core", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("SeuratToScopGiotto", "bridge", "framework_bridge", "GiottoObject.R", "Giotto", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed"),
    entry("RunGiottoWorkflow", "workflow", "framework_workflow", "GiottoObject.R", "Giotto", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoPreprocess", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoReduce", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoCluster", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoSpatialNetwork", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoSpatialGenes", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoSpatialModules", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoCellProximity", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoHMRF", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("AddGiottoToSeurat", "bridge", "framework_bridge", "GiottoObject.R", "Giotto", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed"),
    entry("RunGiottoCluster", "analysis", "framework_workflow", "RunGiottoCluster.R", "GiottoCluster", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("RunGiottoCellProximity", "analysis", "framework_workflow", "RunGiottoSpatialMethods.R", "GiottoCellProximity", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("RunGiottoSpatialGenes", "analysis", "framework_workflow", "RunGiottoSpatialMethods.R", "GiottoSpatialGenes", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("RunGiottoSpatialModules", "analysis", "framework_workflow", "RunGiottoSpatialMethods.R", "GiottoSpatialModules", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),

    entry("RunDeconvolution", "analysis", "deconvolution", "RunDeconvolution.R", backend_id = "music;bisquerna;bayesprism;cibersort", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "DeconvolutionPlot", backend_requirement = "any"),
    entry("RunRCTD", "analysis", "deconvolution", "RunRCTD.R", "RCTD", "spacexr", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialDeconvolutionPlot"),
    entry("RunCSIDE", "analysis", "deconvolution", "RunCSIDE.R", "CSIDE", "spacexr", coordinate_space_current = "mixed", coordinate_space_target = "mixed", coordinate_requirement = "backend_managed"),
    entry("RunCARD", "analysis", "deconvolution", "RunCARD.R", "CARD", "card", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialDeconvolutionPlot"),
    entry("RunSTdeconvolve", "analysis", "deconvolution", "RunSTdeconvolve.R", "STdeconvolve", "stdeconvolve", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "STdeconvolvePlot"),
    entry("RunSPOTlight", "analysis", "deconvolution", "RunSPOTlight.R", "SPOTlight", "spotlight", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "SpatialDeconvolutionPlot"),
    entry("RunCell2location", "analysis", "deconvolution", "RunCell2location.R", "Cell2location", "cell2location", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "Cell2locationPlot"),
    entry("RunSpatialDWLS", "analysis", "deconvolution", "RunSpatialDWLS.R", "SpatialDWLS", "core", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialDeconvolutionPlot"),
    entry("RunSpatialEcoTyper", "analysis", "ecotype", "RunSpatialEcoTyper.R", "SpatialEcoTyper", "spatialecotyper", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "SpatialEcoTyperSpatialPlot"),

    entry("RunBayesSpace", "analysis", "domain", "RunBayesSpace.R", "BayesSpace", "bayesspace", coordinate_space_current = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),
    entry("RunBANKSY", "analysis", "domain", "RunBANKSY.R", "BANKSY", "banksy", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),
    entry("RunCytoSPACE", "analysis", "mapping", "RunCytoSPACE.R", "CytoSPACE", "core", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),
    entry("RunSmoothClust", "analysis", "domain", "RunSmoothClust.R", "SmoothClust", "smoothclust", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", scalability = "sparse_required", plot_function = "SpatialSpotPlot"),
    entry("RunBenchmark", "workflow", "benchmark", "RunBenchmark.R", backend_id = "bayesspace;banksy;smoothclust", coordinate_space_current = "mixed", coordinate_space_target = "mixed", coordinate_requirement = "backend_managed", plot_function = "BenchmarkPlot", backend_requirement = "any"),
    entry("RunMERINGUE", "analysis", "feature_pattern", "RunMERINGUE.R", "MERINGUE", "meringue", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialVariableFeaturePlot"),
    entry("RunSpatialVariableFeatures", "analysis", "feature_pattern", "RunSpatialVariableFeatures.R", "SpatialVariableFeatures", "core;sparkx;nnsvg", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", scalability = "sparse_required", plot_function = "SpatialVariableFeaturePlot", backend_requirement = "any"),
    entry("RunSpatialGradientFeatures", "analysis", "feature_pattern", "RunSpatialGradientFeatures.R", "SpatialGradientFeatures", "core;spata2", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialGradientPlot", backend_requirement = "any"),

    entry("RunSpatialNetwork", "analysis", "network", "RunSpatialNetwork.R", "SpatialNetwork", "biocneighbors", coordinate_space_current = "raw", coordinate_requirement = "distance_sensitive", scalability = "sparse", plot_function = "SpatialNetworkPlot"),
    entry("RunSpatialNeighborhood", "analysis", "neighborhood", "RunSpatialNeighborhood.R", "SpatialNeighborhood", "core;spicyr", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", scalability = "sparse_required", plot_function = "SpatialNeighborhoodPlot", backend_requirement = "any"),
    entry("RunSpatialCellChat", "analysis", "communication", "RunSpatialCellChat.R", "SpatialCellChat", "spatialcellchat", coordinate_space_current = "raw", coordinate_requirement = "distance_sensitive", scalability = "sparse_required", plot_function = "SpatialCellChatPlot"),
    entry("RunStatialKontextual", "analysis", "neighborhood", "RunStatialKontextual.R", "StatialKontextual", "statial", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "StatialKontextualPlot"),
    entry("RunSpatialIntegration", "analysis", "integration", "RunSpatialIntegration.R", "SpatialIntegration", "precast;bass;spatialmnn", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialIntegrationPlot", backend_requirement = "any"),
    entry("RunMistyR", "analysis", "neighborhood", "RunMistyR.R", "MistyR", "mistyr", coordinate_space_current = "raw", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "MistyRPlot"),
    entry("RunSemlaSpatialNetwork", "analysis", "neighborhood", "RunSemla.R", "SemlaSpatialNetwork", "semla", coordinate_space_current = "raw", coordinate_requirement = "distance_sensitive"),
    entry("RunSemlaLocalG", "analysis", "neighborhood", "RunSemla.R", "SemlaLocalG", "semla", coordinate_space_current = "raw", coordinate_requirement = "distance_sensitive"),
    entry("RunSemlaRadialDistance", "analysis", "neighborhood", "RunSemla.R", "SemlaRadialDistance", "semla", coordinate_space_current = "raw", coordinate_requirement = "distance_sensitive"),
    entry("RunSemlaRegionNeighbors", "analysis", "neighborhood", "RunSemla.R", "SemlaRegionNeighbors", "semla", coordinate_space_current = "raw", coordinate_requirement = "distance_sensitive"),

    entry("GiottoPlot", "plot", "visualization", "GiottoPlot.R", backend_id = "giotto;giotto_class", status = "legacy", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("DeconvolutionPlot", "plot", "visualization", "DeconvolutionPlot.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialDeconvolutionPlot", "plot", "visualization", "SpatialDeconvolutionPlot.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialEcoTyperCompositionPlot", "plot", "visualization", "RunSpatialEcoTyper.R", coordinate_space_current = "none"),
    entry("SpatialEcoTyperSpatialPlot", "plot", "visualization", "RunSpatialEcoTyper.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialGradientPlot", "plot", "visualization", "RunSpatialGradientFeatures.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialIntegrationPlot", "plot", "visualization", "RunSpatialIntegration.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialNeighborhoodPlot", "plot", "visualization", "RunSpatialNeighborhood.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialNetworkPlot", "plot", "visualization", "RunSpatialNetwork.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialCellChatPlot", "plot", "visualization", "RunSpatialCellChat.R", backend_id = "core", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialCellPlot", "plot", "visualization", "SpatialCellPlot.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialSpotPlot", "plot", "visualization", "SpatialSpotPlot.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialVariableFeaturePlot", "plot", "visualization", "RunSpatialVariableFeatures.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("MistyRPlot", "plot", "visualization", "RunMistyR.R", coordinate_space_current = "none"),
    entry("StatialKontextualPlot", "plot", "visualization", "RunStatialKontextual.R", coordinate_space_current = "none"),
    entry("STdeconvolvePlot", "plot", "visualization", "RunSTdeconvolve.R", backend_id = "stdeconvolve", coordinate_space_current = "none"),
    entry("Cell2locationPlot", "plot", "visualization", "RunCell2location.R", backend_id = "cell2location", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("BenchmarkPlot", "plot", "visualization", "BenchmarkPlot.R", backend_id = "core", coordinate_space_current = "none"),
    entry("standard_scop", "workflow", "recommended_workflow", "standard_scop.R", backend_id = "core", coordinate_space_current = "mixed", coordinate_space_target = "mixed", coordinate_requirement = "backend_managed")
  )
  registry <- do.call(rbind, rows)
  rownames(registry) <- NULL
  registry
}

spatial_giotto_symbol_registry <- function() {
  class_symbols <- c(
    "createExprObj", "createGiottoObject", "createNearestNetwork",
    "createSpatialNetwork", "fDataDT", "getDimReduction",
    "getSpatialNetwork", "pDataDT", "setExpression"
  )
  giotto_symbols <- c(
    "normalizeGiotto", "runPCA", "runUMAP", "doLeidenCluster",
    "doLouvainCluster", "binSpect", "detectSpatialCorFeats",
    "clusterSpatialCorFeats", "cellProximityEnrichment", "doHMRF",
    "addHMRF", "addStatistics", "calculateHVF"
  )
  data.frame(
    symbol = c(class_symbols, giotto_symbols),
    package = c(
      rep("GiottoClass", length(class_symbols)),
      rep("Giotto", length(giotto_symbols))
    ),
    required = !c(
      rep(FALSE, length(class_symbols)),
      giotto_symbols %in% c("addStatistics", "calculateHVF")
    ),
    stringsAsFactors = FALSE
  )
}

spatial_giotto_converter_name <- function() {
  seurat_major <- suppressWarnings(as.integer(
    strsplit(as.character(utils::packageVersion("Seurat")), "\\.")[[1L]][1L]
  ))
  if (!is.na(seurat_major) && seurat_major >= 5L) {
    "giottoToSeuratV5"
  } else {
    "giottoToSeuratV4"
  }
}

spatial_backend_registry <- function() {
  backend <- function(
    id,
    package,
    repository = package,
    symbols = character(),
    symbol_sets = list(),
    method_symbols = list(),
    runtime = "r",
    environment_modules = character(),
    requirements = character(),
    package_candidates = package
  ) {
    list(
      id = id,
      package = package,
      repository = repository,
      symbols = symbols,
      symbol_sets = symbol_sets,
      method_symbols = method_symbols,
      runtime = runtime,
      environment_modules = environment_modules,
      requirements = requirements,
      package_candidates = unique(package_candidates)
    )
  }
  giotto_symbols <- spatial_giotto_symbol_registry()
  list(
    core = backend("core", "scop", runtime = "core"),
    biocneighbors = backend("biocneighbors", "BiocNeighbors", symbols = c("findKNN", "findNeighbors")),
    spatialcellchat = backend(
      "spatialcellchat",
      "SpatialCellChat",
      "jinworks/SpatialCellChat",
      c(
        "createSpatialCellChat", "subsetDB", "subsetData", "preProcessing",
        "identifyOverExpressedGenes", "identifyOverExpressedInteractions",
        "computeCommunProb", "filterProbability", "filterCommunication",
        "computeAvgCommunProb", "computeAvgCommunProb_Visium",
        "computeCommunProbPathway", "aggregateNet", "subsetCommunication"
      )
    ),
    spanorm = backend("spanorm", "SpaNorm", symbols = "SpaNorm"),
    spatialqm = backend(
      "spatialqm",
      "SpatialQM",
      symbols = c(
        "getNcells", "getTxPerCell", "getTxPerArea", "getTxPerNuc",
        "getMeanExpression", "getMeanSignalRatio", "getCellTxFraction",
        "getMaxRatio", "getMaxDetection", "getMECR", "getMorans",
        "getSilhouetteWidth", "getSparsity", "getEntropy"
      )
    ),
    spotsweeper = backend(
      "spotsweeper",
      "SpotSweeper",
      symbols = c("localOutliers", "findArtifacts")
    ),
    giotto = backend(
      "giotto", "Giotto", "drieslab/Giotto",
      giotto_symbols$symbol[giotto_symbols$package == "Giotto" & giotto_symbols$required]
    ),
    giotto_class = backend(
      "giotto_class", "GiottoClass", "drieslab/GiottoClass",
      giotto_symbols$symbol[giotto_symbols$package == "GiottoClass" & giotto_symbols$required]
    ),
    spata2 = backend("spata2", "SPATA2", "theMILOlab/SPATA2", c("asSPATA2", "asSeurat")),
    music = backend("music", "MuSiC", "xuranw/MuSiC", "music_prop"),
    bisquerna = backend("bisquerna", "BisqueRNA", symbols = "ReferenceBasedDecomposition"),
    bayesprism = backend("bayesprism", "BayesPrism", symbols = c("new.prism", "run.prism")),
    cibersort = backend("cibersort", "CIBERSORT", "Moonerss/CIBERSORT", "cibersort"),
    spacexr = backend(
      "spacexr",
      "spacexr",
      "dmcable/spacexr",
      symbol_sets = list(
        new = c("createRctd", "runRctd"),
        legacy = c("SpatialRNA", "Reference", "create.RCTD", "run.RCTD")
      )
    ),
    card = backend(
      "card",
      "CARD",
      "YingMa0107/CARD",
      c("createCARDObject", "CARD_deconvolution"),
      package_candidates = c("CARD", "CARDspa")
    ),
    stdeconvolve = backend("stdeconvolve", "STdeconvolve", "JEFworks-Lab/STdeconvolve", c("cleanCounts", "fitLDA")),
    spotlight = backend("spotlight", "SPOTlight", symbols = "SPOTlight"),
    cell2location = backend(
      "cell2location",
      "cell2location",
      "BayraktarLab/cell2location",
      runtime = "python",
      environment_modules = c("cell2location", "scanpy", "scvi"),
      requirements = c(
        "cell2location==0.1.5", "scvi-tools==1.3.3", "scanpy",
        "anndata", "numpy", "pandas", "scipy", "torch"
      )
    ),
    spatialecotyper = backend(
      "spatialecotyper",
      "SpatialEcoTyper",
      "digitalcytometry/SpatialEcoTyper",
      c("SpatialEcoTyper", "MultiSpatialEcoTyper", "RecoverSE", "DeconvoluteSE")
    ),
    bayesspace = backend(
      "bayesspace",
      "BayesSpace",
      symbols = c("spatialPreprocess", "spatialCluster")
    ),
    banksy = backend("banksy", "Banksy", symbols = "computeBanksy"),
    smoothclust = backend(
      "smoothclust",
      "smoothclust",
      "lmweber/smoothclust",
      "smoothclust"
    ),
    meringue = backend(
      "meringue",
      "MERINGUE",
      symbols = c(
        "getSpatialNeighbors", "moranTest", "moranPermutationTest",
        "spatialCrossCorMatrix", "spatialCrossCorTest",
        "groupSigSpatialPatterns"
      )
    ),
    sparkx = backend("sparkx", "SPARK", symbols = "sparkx"),
    nnsvg = backend("nnsvg", "nnSVG", symbols = "nnSVG"),
    spicyr = backend("spicyr", "spicyR", symbols = "spicy"),
    statial = backend("statial", "Statial", symbols = "Kontextual"),
    precast = backend(
      "precast",
      "PRECAST",
      "feiyoung/PRECAST",
      c("CreatePRECASTObject", "AddAdjList", "AddParSetting", "PRECAST", "SelectModel")
    ),
    bass = backend(
      "bass",
      "BASS",
      "zhengli09/BASS",
      symbols = c("createBASSObject", "BASS.preprocess", "BASS.run")
    ),
    spatialmnn = backend(
      "spatialmnn",
      "atlasClustering",
      "Pixel-Dream/spatialMNN",
      c("stage_1", "stage_2")
    ),
    mistyr = backend("mistyr", "mistyR", symbols = c("create_initial_view", "run_misty")),
    semla = backend(
      "semla",
      "semla",
      symbols = "UpdateSeuratForSemla",
      method_symbols = list(
        RunSemlaSpatialNetwork = "GetSpatialNetwork",
        RunSemlaLocalG = "RunLocalG",
        RunSemlaRegionNeighbors = "RegionNeighbors",
        RunSemlaRadialDistance = "RadialDistance"
      )
    )
  )
}

spatial_registry_split_backends <- function(x) {
  unique(trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)))
}

#' @title List spatial methods
#'
#' @description
#' List the public spatial analysis, visualization, framework bridge, and
#' workflow entry points registered by scop. This function does not inspect
#' optional packages unless `available` is `TRUE` or `FALSE`.
#'
#' @param task,kind,backend,status Optional exact filters.
#' @param available Optional logical availability filter. `NULL` avoids backend
#'   inspection.
#' @param pattern Optional case-insensitive pattern matched against method,
#'   task, backend, and documentation fields.
#'
#' @return A data frame with one row per registered public method.
#'
#' @export
ListSpatialMethods <- function(
  task = NULL,
  kind = NULL,
  backend = NULL,
  status = NULL,
  available = NULL,
  pattern = NULL
) {
  registry <- spatial_method_registry()
  exact_filter <- function(data, column, value) {
    if (is.null(value)) return(data)
    data[data[[column]] %in% value, , drop = FALSE]
  }
  registry <- exact_filter(registry, "task", task)
  registry <- exact_filter(registry, "kind", kind)
  registry <- exact_filter(registry, "status", status)
  if (!is.null(backend)) {
    keep <- vapply(
      registry$backend_id,
      function(x) any(spatial_registry_split_backends(x) %in% backend),
      logical(1)
    )
    registry <- registry[keep, , drop = FALSE]
  }
  if (!is.null(pattern)) {
    if (!is.character(pattern) || length(pattern) != 1L || is.na(pattern) || !nzchar(pattern)) {
      log_message("{.arg pattern} must be one non-empty string", message_type = "error")
    }
    haystack <- apply(
      registry[, c("method", "task", "backend_id", "documentation"), drop = FALSE],
      1,
      paste,
      collapse = " "
    )
    registry <- registry[grepl(pattern, haystack, ignore.case = TRUE), , drop = FALSE]
  }
  if (!is.null(available)) {
    if (!is.logical(available) || length(available) != 1L || is.na(available)) {
      log_message("{.arg available} must be TRUE, FALSE, or NULL", message_type = "error")
    }
    status_table <- SpatialBackendStatus(api_check = TRUE)
    ready <- stats::setNames(status_table$availability == "available", status_table$backend_id)
    method_ready <- vapply(seq_len(nrow(registry)), function(i) {
      ids <- spatial_registry_split_backends(registry$backend_id[[i]])
      values <- ready[ids]
      values[is.na(values)] <- FALSE
      if (identical(registry$backend_requirement[[i]], "any")) {
        any(values)
      } else {
        all(values)
      }
    }, logical(1))
    registry <- registry[method_ready == available, , drop = FALSE]
  }
  rownames(registry) <- NULL
  registry
}

spatial_backend_cache_key <- function(api_check, backend_ids, method, envname, conda) {
  paste(
    "v1",
    paste(normalizePath(.libPaths(), winslash = "/", mustWork = FALSE), collapse = "|"),
    isTRUE(api_check),
    paste(sort(backend_ids), collapse = "|"),
    paste(sort(method %||% character()), collapse = "|"),
    envname,
    conda,
    sep = "::"
  )
}

spatial_backend_readonly_conda <- function(conda = "auto") {
  if (!is.character(conda) || length(conda) != 1L || is.na(conda) || !nzchar(conda)) {
    log_message("{.arg conda} must be one non-empty string", message_type = "error")
  }
  if (identical(conda, "auto")) {
    return(tryCatch(find_conda(), error = function(e) NULL))
  }
  tryCatch(
    resolve_conda_executable(
      conda,
      error_if_missing = FALSE,
      install_if_missing = FALSE
    ),
    error = function(e) NULL
  )
}

spatial_python_backend_status <- function(spec, envname, conda, api_check) {
  requirements <- spec$requirements %||% character()
  conda_path <- spatial_backend_readonly_conda(conda)
  namespace_error <- NA_character_
  environment_exists <- FALSE
  if (!is.null(conda_path)) {
    environment_exists <- tryCatch(
      isTRUE(env_exist(conda = conda_path, envname = envname)),
      error = function(e) {
        namespace_error <<- conditionMessage(e)
        FALSE
      }
    )
  }
  package_status <- stats::setNames(rep(FALSE, length(requirements)), requirements)
  if (environment_exists && isTRUE(api_check) && length(requirements) > 0L) {
    package_status <- tryCatch(
      exist_python_pkgs(
        packages = requirements,
        envname = envname,
        conda = conda_path,
        verbose = FALSE
      ),
      error = function(e) {
        namespace_error <<- conditionMessage(e)
        stats::setNames(rep(FALSE, length(requirements)), requirements)
      }
    )
  }
  missing <- if (isTRUE(api_check)) names(package_status)[!package_status] else character()
  installed <- environment_exists && (!isTRUE(api_check) || length(missing) == 0L)
  availability <- if (!environment_exists) {
    "missing"
  } else if (!is.na(namespace_error)) {
    "namespace_error"
  } else if (length(missing) > 0L) {
    "api_incompatible"
  } else {
    "available"
  }
  data.frame(
    backend_id = spec$id,
    runtime = spec$runtime,
    environment = envname,
    package = spec$package,
    repository = spec$repository,
    installed = installed,
    api_checked = isTRUE(api_check) && length(requirements) > 0L,
    required_symbols = paste(requirements, collapse = ","),
    missing_symbols = paste(missing, collapse = ","),
    availability = availability,
    namespace_error = namespace_error,
    stringsAsFactors = FALSE
  )
}

spatial_backend_required_symbols <- function(spec, method = NULL, exports = NULL) {
  symbols <- spec$symbols
  symbol_sets <- spec$symbol_sets %||% list()
  method_symbols <- spec$method_symbols %||% list()
  if (length(symbol_sets) > 0L) {
    if (is.null(exports)) {
      symbols <- symbol_sets[[1L]]
    } else {
      complete <- vapply(
        symbol_sets,
        function(candidate) all(candidate %in% exports),
        logical(1)
      )
      if (any(complete)) {
        symbols <- symbol_sets[[which(complete)[[1L]]]]
      } else {
        matched <- vapply(
          symbol_sets,
          function(candidate) sum(candidate %in% exports),
          integer(1)
        )
        symbols <- symbol_sets[[which.max(matched)]]
      }
    }
  }
  if (length(method_symbols) > 0L) {
    selected <- if (is.null(method)) {
      unlist(method_symbols, use.names = FALSE)
    } else {
      matched <- intersect(as.character(method), names(method_symbols))
      unlist(method_symbols[matched], use.names = FALSE)
    }
    symbols <- c(symbols, selected)
  }
  if (
    identical(spec$id, "giotto_class") &&
      !is.null(method) &&
      "giotto_to_srt" %in% method
  ) {
    symbols <- c(symbols, spatial_giotto_converter_name())
  }
  unique(symbols)
}

spatial_backend_resolve_package <- function(spec, installed = NULL) {
  candidates <- unique(spec$package_candidates %||% spec$package)
  if (is.null(installed)) {
    installed <- rownames(utils::installed.packages())
  }
  available <- candidates[candidates %in% installed]
  if (length(available) > 0L) available[[1L]] else spec$package
}

#' @title Inspect spatial backend availability
#'
#' @description
#' Diagnose optional spatial backend installation and exported API status. The
#' function is read-only: it never installs or updates a package. Actual method
#' calls remain the authoritative runtime check for system libraries, S4
#' registration, and data-specific requirements.
#'
#' @param backend Optional backend identifier filter.
#' @param method Optional registered spatial method filter.
#' @param api_check Whether to inspect required namespace exports.
#' @param refresh Whether to refresh the session-local diagnostic cache.
#' @param envname Existing SCOP Python environment to inspect for Python
#'   backends. The diagnostic never creates or modifies the environment.
#' @param conda Existing conda-compatible executable, or `"auto"` to discover
#'   one without installing it.
#'
#' @return A data frame with one row per backend.
#'
#' @export
SpatialBackendStatus <- function(
  backend = NULL,
  method = NULL,
  api_check = TRUE,
  refresh = FALSE,
  envname = NULL,
  conda = "auto"
) {
  if (!is.logical(api_check) || length(api_check) != 1L || is.na(api_check)) {
    log_message("{.arg api_check} must be TRUE or FALSE", message_type = "error")
  }
  if (!is.logical(refresh) || length(refresh) != 1L || is.na(refresh)) {
    log_message("{.arg refresh} must be TRUE or FALSE", message_type = "error")
  }
  backends <- spatial_backend_registry()
  envname <- get_envname(envname)
  backend_ids <- names(backends)
  if (!is.null(method)) {
    registry <- spatial_method_registry()
    registry <- registry[registry$method %in% method, , drop = FALSE]
    backend_ids <- intersect(
      backend_ids,
      spatial_registry_split_backends(registry$backend_id)
    )
  }
  if (!is.null(backend)) {
    backend_ids <- intersect(backend_ids, backend)
  }
  backends <- backends[backend_ids]
  if (length(backends) == 0L) {
    return(data.frame(
      backend_id = character(), runtime = character(), environment = character(),
      package = character(), repository = character(),
      installed = logical(), api_checked = logical(), required_symbols = character(),
      missing_symbols = character(), availability = character(), namespace_error = character(),
      stringsAsFactors = FALSE
    ))
  }
  cache_key <- spatial_backend_cache_key(api_check, backend_ids, method, envname, conda)
  if (isTRUE(refresh) && exists(cache_key, envir = .spatial_backend_cache, inherits = FALSE)) {
    rm(list = cache_key, envir = .spatial_backend_cache)
  }
  if (!exists(cache_key, envir = .spatial_backend_cache, inherits = FALSE)) {
    installed <- rownames(utils::installed.packages())
    rows <- lapply(backends, function(spec) {
      if (identical(spec$runtime, "python")) {
        return(spatial_python_backend_status(
          spec = spec,
          envname = envname,
          conda = conda,
          api_check = api_check
        ))
      }
      package <- spatial_backend_resolve_package(spec, installed = installed)
      symbols <- spatial_backend_required_symbols(spec, method = method)
      is_core <- identical(spec$runtime, "core")
      is_installed <- is_core || package %in% installed
      namespace_error <- NA_character_
      missing_symbols <- character()
      if (!is_installed) {
        availability <- "missing"
      } else if (isTRUE(api_check) && length(symbols) > 0L) {
        namespace_symbols <- tryCatch(
          ls(asNamespace(package), all.names = TRUE),
          error = function(e) {
            namespace_error <<- conditionMessage(e)
            NULL
          }
        )
        if (is.null(namespace_symbols)) {
          availability <- "namespace_error"
        } else {
          symbols <- spatial_backend_required_symbols(
            spec,
            method = method,
            exports = namespace_symbols
          )
          missing_symbols <- setdiff(symbols, namespace_symbols)
          availability <- if (length(missing_symbols) == 0L) "available" else "api_incompatible"
        }
      } else {
        availability <- "available"
      }
      data.frame(
        backend_id = spec$id,
        runtime = spec$runtime,
        environment = NA_character_,
        package = package,
        repository = spec$repository,
        installed = is_installed,
        api_checked = isTRUE(api_check) && length(symbols) > 0L,
        required_symbols = paste(symbols, collapse = ","),
        missing_symbols = paste(missing_symbols, collapse = ","),
        availability = availability,
        namespace_error = namespace_error,
        stringsAsFactors = FALSE
      )
    })
    assign(cache_key, do.call(rbind, rows), envir = .spatial_backend_cache)
  }
  out <- get(cache_key, envir = .spatial_backend_cache, inherits = FALSE)
  rownames(out) <- NULL
  out
}

spatial_result_backend_versions <- function(backend_id) {
  ids <- spatial_registry_split_backends(backend_id %||% character())
  registry <- spatial_backend_registry()
  selected <- registry[intersect(ids, names(registry))]
  selected <- selected[vapply(selected, function(spec) identical(spec$runtime, "r"), logical(1))]
  installed <- rownames(utils::installed.packages())
  packages <- vapply(
    selected,
    spatial_backend_resolve_package,
    character(1),
    installed = installed
  )
  versions <- vapply(packages, function(package) {
    if (identical(package, "scop")) return(NA_character_)
    tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
  }, character(1))
  versions[!is.na(versions)]
}

spatial_result_build <- function(
  bundle = list(),
  method,
  result_type,
  source = list(),
  provenance = list(),
  parameters = NULL,
  summary = NULL
) {
  if (!is.list(bundle)) bundle <- list(result = bundle)
  if (!is.character(method) || length(method) != 1L || is.na(method) || !nzchar(method)) {
    log_message("{.arg method} must be one non-empty result family", message_type = "error")
  }
  if (!is.character(result_type) || length(result_type) != 1L || is.na(result_type) || !nzchar(result_type)) {
    log_message("{.arg result_type} must be one non-empty registry task", message_type = "error")
  }
  if (!is.list(source) || !is.list(provenance)) {
    log_message("{.arg source} and {.arg provenance} must be lists", message_type = "error")
  }
  source_input <- source
  parameters <- parameters %||% bundle$parameters %||% list()
  summary <- summary %||% bundle$summary %||% list()
  if (!is.list(parameters)) {
    log_message("Spatial result {.arg parameters} must be a list", message_type = "error")
  }
  coordinate_payload <- bundle$coords %||% bundle$spatial_coords %||% NULL
  coordinate_source <- attr(coordinate_payload, "spatial_source") %||% list()
  coordinate_transform <- attr(coordinate_payload, "spatial_transform")
  source_defaults <- list(
    image = parameters$image %||% character(),
    coordinate_space = parameters$coordinate_space %||% "none",
    transform = coordinate_transform,
    assay = parameters$assay %||% NA_character_,
    layer = parameters$layer %||% NA_character_
  )
  source <- utils::modifyList(source_defaults, coordinate_source)
  source <- utils::modifyList(source, source_input)
  provenance_defaults <- list(
    producer = NA_character_, backend_id = "core",
    backend_versions = character(),
    scop_version = tryCatch(as.character(utils::packageVersion("scop")), error = function(e) NA_character_)
  )
  provenance <- utils::modifyList(provenance_defaults, provenance)
  if (length(provenance$backend_versions) == 0L) {
    provenance$backend_versions <- spatial_result_backend_versions(provenance$backend_id)
  }
  bundle$method <- method
  bundle$schema_version <- 1L
  bundle$result_type <- result_type
  bundle$source <- source
  bundle$provenance <- provenance
  bundle$parameters <- parameters
  bundle$summary <- summary
  spatial_result_validate(bundle)
  bundle
}

spatial_result_validate <- function(bundle) {
  required <- c(
    "method", "schema_version", "result_type", "source",
    "provenance", "parameters", "summary"
  )
  if (!is.list(bundle) || !all(required %in% names(bundle))) {
    log_message(
      "Spatial result schema v1 requires fields: {.val {required}}",
      message_type = "error"
    )
  }
  if (!identical(as.integer(bundle$schema_version), 1L)) {
    log_message("Unsupported spatial result schema version", message_type = "error")
  }
  if (!is.list(bundle$source) || !is.list(bundle$provenance) || !is.list(bundle$parameters)) {
    log_message("Spatial result source, provenance, and parameters must be lists", message_type = "error")
  }
  invisible(TRUE)
}

spatial_result_registry_row <- function(tool_name = NULL, bundle = NULL) {
  registry <- spatial_method_registry()
  registry <- registry[!is.na(registry$tool_key) & nzchar(registry$tool_key), , drop = FALSE]
  producer <- if (is.list(bundle)) bundle$provenance$producer %||% NA_character_ else NA_character_
  family <- if (is.list(bundle)) bundle$method %||% NA_character_ else NA_character_
  match_index <- integer()
  if (!is.null(tool_name)) match_index <- which(registry$tool_key == tool_name)
  if (length(match_index) == 0L && length(producer) == 1L && !is.na(producer)) {
    match_index <- which(registry$method == producer)
  }
  if (length(match_index) == 0L && length(family) == 1L && !is.na(family)) {
    match_index <- which(registry$tool_key == family)
  }
  if (length(match_index) == 0L) return(NULL)
  registry[match_index[[1L]], , drop = FALSE]
}

spatial_result_normalize <- function(bundle, tool_name = NULL, registry_row = NULL) {
  if (is.list(bundle) && identical(as.integer(bundle$schema_version %||% NA_integer_), 1L)) {
    spatial_result_validate(bundle)
    return(bundle)
  }
  registry_row <- registry_row %||% spatial_result_registry_row(tool_name, bundle)
  if (is.null(registry_row)) {
    log_message("Stored tool {.val {tool_name}} is not a registered spatial result", message_type = "error")
  }
  parameters <- if (is.list(bundle)) bundle$parameters %||% list() else list()
  source <- if (is.list(bundle)) bundle$source %||% list() else list()
  source <- utils::modifyList(
    list(
      image = parameters$image %||% character(),
      coordinate_space = parameters$coordinate_space %||% registry_row$coordinate_space_current[[1L]],
      transform = NULL,
      assay = parameters$assay %||% NA_character_,
      layer = parameters$layer %||% NA_character_
    ),
    source
  )
  spatial_result_build(
    bundle = bundle,
    method = registry_row$tool_key[[1L]],
    result_type = registry_row$task[[1L]],
    source = source,
    provenance = list(
      producer = registry_row$method[[1L]],
      backend_id = registry_row$backend_id[[1L]]
    ),
    parameters = parameters,
    summary = if (is.list(bundle)) bundle$summary %||% list() else list()
  )
}

spatial_result_index <- function(object) {
  registry <- spatial_method_registry()
  registry <- registry[!is.na(registry$tool_key) & nzchar(registry$tool_key), , drop = FALSE]
  keys <- names(object@tools)
  rows <- lapply(keys, function(key) {
    bundle <- object@tools[[key]]
    row <- spatial_result_registry_row(key, bundle)
    if (is.null(row)) return(NULL)
    is_schema <- is.list(bundle) && identical(as.integer(bundle$schema_version %||% NA_integer_), 1L)
    is_default <- key %in% registry$tool_key
    if (!is_schema && !is_default) return(NULL)
    data.frame(
      tool_name = key,
      registry_method = row$method[[1L]],
      family = row$tool_key[[1L]],
      result_type = row$task[[1L]],
      backend_id = row$backend_id[[1L]],
      plot_function = row$plot_function[[1L]],
      coordinate_space_current = row$coordinate_space_current[[1L]],
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      tool_name = character(), registry_method = character(), family = character(),
      result_type = character(), backend_id = character(), plot_function = character(),
      coordinate_space_current = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

spatial_result_payload_size <- function(x) {
  if (is.null(x)) return(0L)
  if (is.data.frame(x) || is.matrix(x)) return(nrow(x))
  if (is.atomic(x)) return(length(x))
  if (is.list(x)) return(length(x))
  1L
}

spatial_result_state <- function(bundle, object_cells = NULL) {
  metadata_names <- c(
    "method", "schema_version", "result_type", "source", "provenance",
    "parameters", "summary", "active_graph", "active_method"
  )
  if (!is.list(bundle)) {
    return(list(state = if (length(bundle) > 0L) "ready" else "empty", n_items = length(bundle), reason = NA_character_))
  }
  if (!is.null(bundle$graphs) && length(bundle$graphs) > 0L) {
    node_ids <- unique(unlist(lapply(bundle$graphs, function(graph) {
      nodes <- graph$nodes
      if (is.data.frame(nodes) && "cell_id" %in% colnames(nodes)) as.character(nodes$cell_id) else character()
    }), use.names = FALSE))
    stale <- length(node_ids) > 0L && !is.null(object_cells) && any(!node_ids %in% object_cells)
    return(list(
      state = if (stale) "stale" else "ready",
      n_items = length(bundle$graphs),
      reason = if (stale) "stored graph nodes are absent from object" else NA_character_
    ))
  }
  payload_ids <- character()
  if (is.data.frame(bundle$nodes) && "cell_id" %in% colnames(bundle$nodes)) {
    payload_ids <- c(payload_ids, as.character(bundle$nodes$cell_id))
  }
  if (is.atomic(bundle$cells)) {
    payload_ids <- c(payload_ids, as.character(bundle$cells))
  }
  payload_ids <- unique(payload_ids[!is.na(payload_ids) & nzchar(payload_ids)])
  if (length(payload_ids) > 0L && !is.null(object_cells) && any(!payload_ids %in% object_cells)) {
    return(list(
      state = "stale",
      n_items = length(payload_ids),
      reason = "stored cells or nodes are absent from object"
    ))
  }
  payload_names <- setdiff(names(bundle), metadata_names)
  payload_sizes <- vapply(bundle[payload_names], spatial_result_payload_size, integer(1))
  n_items <- sum(payload_sizes)
  if (n_items == 0L) {
    return(list(state = "empty", n_items = 0L, reason = "no logical result payload"))
  }
  partial <- is.null(bundle$method) || is.null(bundle$source) || is.null(bundle$parameters)
  list(
    state = if (partial) "partial" else "ready",
    n_items = n_items,
    reason = if (partial) "result provenance is incomplete" else NA_character_
  )
}

#' @title Read one stored spatial result
#'
#' @description
#' Return a schema-v1 read-only view of a registered spatial result. Legacy
#' results are normalized in the returned copy and are never written back.
#'
#' @param object A `Seurat` object.
#' @param method Optional public producer or result family.
#' @param tool_name Optional exact key in `object@tools`.
#' @param raw Whether to return the stored value without normalization.
#' @param validate Whether to validate schema-v1 results before returning.
#'
#' @return A plain spatial result list.
#' @export
GetSpatialResult <- function(
  object,
  method = NULL,
  tool_name = NULL,
  raw = FALSE,
  validate = TRUE
) {
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  }
  if (is.null(method) && is.null(tool_name)) {
    log_message("Provide {.arg method} or {.arg tool_name}", message_type = "error")
  }
  index <- spatial_result_index(object)
  keep <- rep(TRUE, nrow(index))
  if (!is.null(method)) {
    keep <- keep & (
      index$registry_method %in% method | index$family %in% method
    )
  }
  if (!is.null(tool_name)) keep <- keep & index$tool_name %in% tool_name
  matches <- index[keep, , drop = FALSE]
  if (nrow(matches) == 0L) {
    log_message("No stored spatial result matches the requested selector", message_type = "error")
  }
  if (nrow(matches) > 1L) {
    log_message(
      "Multiple spatial results match; select one with {.arg tool_name}: {.val {matches$tool_name}}",
      message_type = "error"
    )
  }
  key <- matches$tool_name[[1L]]
  bundle <- object@tools[[key]]
  if (isTRUE(raw)) return(bundle)
  normalized <- spatial_result_normalize(bundle, tool_name = key)
  if (isTRUE(validate)) spatial_result_validate(normalized)
  normalized
}

#' @title Inspect stored spatial results
#'
#' @description
#' Read spatial result bundles from a Seurat object's `tools` slot without
#' migrating, renaming, or modifying legacy results.
#'
#' @param object A `Seurat` object.
#' @param method Optional registered producer or result-family filter.
#' @param tool_name Optional stored tool key filter.
#' @param include_empty Whether to include recognized keys without a logical
#'   result payload.
#' @param detail Return one row per result or one row per stored spatial graph.
#'
#' @return A data frame describing recognized stored spatial results or graphs.
#' @export
SpatialResultInfo <- function(
  object,
  method = NULL,
  tool_name = NULL,
  include_empty = FALSE,
  detail = c("results", "graphs")
) {
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.logical(include_empty) || length(include_empty) != 1L || is.na(include_empty)) {
    log_message("{.arg include_empty} must be TRUE or FALSE", message_type = "error")
  }
  detail <- match.arg(detail)
  index <- spatial_result_index(object)
  if (!is.null(method)) {
    index <- index[index$registry_method %in% method | index$family %in% method, , drop = FALSE]
  }
  if (!is.null(tool_name)) index <- index[index$tool_name %in% tool_name, , drop = FALSE]
  cells <- colnames(object)
  if (identical(detail, "graphs")) {
    rows <- lapply(seq_len(nrow(index)), function(i) {
      key <- index$tool_name[[i]]
      bundle <- object@tools[[key]]
      if (!is.list(bundle) || is.null(bundle$graphs)) return(NULL)
      graph_names <- names(bundle$graphs)
      if (is.null(graph_names)) graph_names <- as.character(seq_along(bundle$graphs))
      do.call(rbind, lapply(seq_along(bundle$graphs), function(j) {
        graph <- bundle$graphs[[j]]
        source <- graph$source %||% list()
        parameters <- graph$parameters %||% list()
        data.frame(
          tool_name = key,
          graph_name = graph_names[[j]],
          active = identical(graph_names[[j]], bundle$active_graph),
          n_nodes = if (is.data.frame(graph$nodes)) nrow(graph$nodes) else 0L,
          n_edges = if (is.data.frame(graph$edges)) nrow(graph$edges) else 0L,
          method = parameters$method %||% NA_character_,
          weight = parameters$weight %||% NA_character_,
          image = paste(source$image %||% character(), collapse = ","),
          coordinate_space = source$coordinate_space %||% NA_character_,
          stringsAsFactors = FALSE
        )
      }))
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0L) {
      return(data.frame(
        tool_name = character(), graph_name = character(), active = logical(),
        n_nodes = integer(), n_edges = integer(), method = character(), weight = character(),
        image = character(), coordinate_space = character(), stringsAsFactors = FALSE
      ))
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }
  rows <- lapply(seq_len(nrow(index)), function(i) {
    key <- index$tool_name[[i]]
    stored <- object@tools[[key]]
    state <- spatial_result_state(stored, object_cells = cells)
    if (!include_empty && identical(state$state, "empty")) return(NULL)
    normalized <- spatial_result_normalize(stored, tool_name = key)
    source <- normalized$source
    data.frame(
      method = normalized$method,
      producer = normalized$provenance$producer %||% index$registry_method[[i]],
      result_type = normalized$result_type,
      tool_name = key,
      schema_version = if (is.list(stored)) as.integer(stored$schema_version %||% NA_integer_) else NA_integer_,
      result_state = state$state,
      n_items = state$n_items,
      coordinate_space = as.character(source$coordinate_space %||% index$coordinate_space_current[[i]])[[1L]],
      image = paste(source$image %||% character(), collapse = ","),
      parameters_available = is.list(stored) && !is.null(stored$parameters),
      summary_available = is.list(stored) && !is.null(stored$summary),
      plot_function = index$plot_function[[i]],
      empty_reason = state$reason,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      method = character(), producer = character(), result_type = character(),
      tool_name = character(), schema_version = integer(), result_state = character(),
      n_items = integer(), coordinate_space = character(), image = character(),
      parameters_available = logical(), summary_available = logical(),
      plot_function = character(), empty_reason = character(), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
