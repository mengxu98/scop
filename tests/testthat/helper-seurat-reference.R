seurat_reference_method <- function(generic, class, object, ...) {
  method <- get(paste0(generic, ".", class), asNamespace("Seurat"))
  method(object = object, ...)
}

seurat_reference_find_all_markers <- local({
  method <- get("FindAllMarkers", asNamespace("Seurat"))
  method_env <- new.env(parent = environment(method))
  method_env$FindMarkers <- function(object, ...) {
    # Seurat 5.6 also calls FindMarkers on assays with explicit cell groups.
    # Follow S3 inheritance (including Assay5 -> StdAssay) while keeping the
    # reference independent of scop's registered marker methods.
    for (class in c(.class2(object), "default")) {
      marker_method <- get0(
        paste0("FindMarkers.", class),
        envir = asNamespace("Seurat"),
        mode = "function",
        inherits = FALSE
      )
      if (!is.null(marker_method)) {
        return(marker_method(object = object, ...))
      }
    }
    stop("No Seurat FindMarkers method for the reference object.")
  }
  environment(method) <- method_env
  method
})

seurat_reference_add_module_score <- function(object, ...) {
  get("AddModuleScore.Seurat", asNamespace("Seurat"))(
    object = object,
    ...
  )
}

seurat_reference_cell_cycle_scoring <- local({
  method <- get("CellCycleScoring", asNamespace("Seurat"))
  method_env <- new.env(parent = environment(method))
  method_env$AddModuleScore <- seurat_reference_add_module_score
  environment(method) <- method_env
  method
})
