create_xenium_human_pancreas_sub_data <- function(
  datasets_rds = file.path("..", "datasets", "Xenium", "xenium_human_pancreas_sub.rds")
) {
  if (!requireNamespace("usethis", quietly = TRUE)) {
    stop("Package 'usethis' is required to write package data.", call. = FALSE)
  }
  if (!file.exists(datasets_rds)) {
    stop(
      "Missing Xenium source RDS: ",
      normalizePath(datasets_rds, mustWork = FALSE),
      call. = FALSE
    )
  }
  xenium_human_pancreas_sub <- readRDS(datasets_rds)
  if (!inherits(xenium_human_pancreas_sub, "Seurat")) {
    stop("xenium_human_pancreas_sub must be a Seurat object.", call. = FALSE)
  }
  required_meta <- c(
    "x",
    "y",
    "nCount_Xenium",
    "nFeature_Xenium",
    "platform",
    "source_dataset"
  )
  missing_meta <- setdiff(required_meta, colnames(xenium_human_pancreas_sub@meta.data))
  if (length(missing_meta) > 0) {
    stop(
      "Missing required Xenium metadata columns: ",
      paste(missing_meta, collapse = ", "),
      call. = FALSE
    )
  }
  if (!"Xenium" %in% names(xenium_human_pancreas_sub@assays)) {
    stop("xenium_human_pancreas_sub must contain a Xenium assay.", call. = FALSE)
  }
  if (!"TENxXeniumData" %in% names(xenium_human_pancreas_sub@tools)) {
    stop("xenium_human_pancreas_sub must record TENxXeniumData provenance in @tools.", call. = FALSE)
  }
  usethis::use_data(
    xenium_human_pancreas_sub,
    compress = "xz",
    overwrite = TRUE
  )
  invisible(xenium_human_pancreas_sub)
}
