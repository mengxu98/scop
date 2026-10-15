atac_scdblfinder <- function(object, assay, sample.by, db_rate, seed, verbose, ...) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message("{paste(report$errors, collapse = ' ')}", message_type = "error")
  cells <- colnames(object[[assay]])
  metadata <- object[[]][cells, , drop = FALSE]
  if (!is.character(sample.by) || length(sample.by) != 1L ||
    is.na(sample.by) || !sample.by %in% colnames(metadata)) {
    log_message("ATAC doublet calling requires sample.by naming the capture/library metadata column.", message_type = "error")
  }
  captures <- as.character(metadata[[sample.by]])
  if (anyNA(captures) || any(!nzchar(captures))) {
    log_message("Every selected ATAC cell must have a non-empty capture/library ID.", message_type = "error")
  }
  if (!is.null(db_rate) && (!is.numeric(db_rate) || length(db_rate) != 1L ||
    !is.finite(db_rate) || db_rate <= 0 || db_rate >= 1)) {
    log_message("db_rate must be NULL or one probability strictly between zero and one.", message_type = "error")
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
    seed < 0 || seed > .Machine$integer.max || seed %% 1 != 0) {
    log_message("seed must be one nonnegative integer within the R integer range.", message_type = "error")
  }
  dots <- list(...)
  reserved <- c("sce", "samples", "dbr", "aggregateFeatures", "processing", "multiSampleMode", "returnType", "verbose")
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
    anyDuplicated(names(dots)) || any(names(dots) %in% reserved))) {
    log_message("Use named native model arguments; capture grouping, aggregation, processing and return type are fixed by the ATAC contract.", message_type = "error")
  }
  fit <- atac_backend_fun("scDblFinder", "scDblFinder")
  make_sce <- atac_backend_fun("SingleCellExperiment", "SingleCellExperiment")
  with_seed <- atac_backend_fun("withr", "with_seed")
  counts <- SeuratObject::LayerData(object[[assay]], layer = "counts")
  sce <- make_sce(assays = list(counts = counts))
  if (is.null(dots$nfeatures)) dots$nfeatures <- 25L
  if (is.null(dots$BPPARAM)) {
    dots$BPPARAM <- atac_backend_fun("BiocParallel", "SerialParam")(
      RNGseed = seed, progressbar = verbose
    )
  }
  parameters <- c(list(
    samples = captures, dbr = db_rate,
    aggregateFeatures = TRUE, processing = "normFeatures",
    multiSampleMode = "split", returnType = "sce", verbose = verbose
  ), dots)
  fitted <- with_seed(seed, do.call(fit, c(list(sce = sce), parameters)))
  if (!identical(colnames(fitted), cells)) {
    log_message("Native scDblFinder changed the selected cell order.", message_type = "error")
  }
  native <- as.data.frame(SummarizedExperiment::colData(fitted))
  for (name in c("score", "class")) {
    column <- paste0("db.scDblFinder_", name)
    values <- if (column %in% colnames(object[[]])) object[[]][[column]] else rep(NA, ncol(object))
    if (name == "class") values <- as.character(values)
    values[match(cells, colnames(object))] <- if (name == "class") {
      as.character(native[[paste0("scDblFinder.", name)]])
    } else {
      native[[paste0("scDblFinder.", name)]]
    }
    object[[column]] <- values
  }
  SeuratObject::Misc(object[[assay]], "scop_atac_doublets") <- list(
    cells = cells, captures = stats::setNames(captures, cells), sample.by = sample.by,
    parameters = parameters, seed = seed, native = native,
    backend = c(scDblFinder = as.character(utils::packageVersion("scDblFinder")))
  )
  object
}
