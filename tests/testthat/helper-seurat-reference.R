seurat_reference_method <- function(generic, class, object, ...) {
  method <- get(paste0(generic, ".", class), asNamespace("Seurat"))
  method(object = object, ...)
}

seurat_reference_find_all_markers <- local({
  method <- get("FindAllMarkers", asNamespace("Seurat"))
  method_env <- new.env(parent = environment(method))
  method_env$FindMarkers <- function(object, ...) {
    get("FindMarkers.Seurat", asNamespace("Seurat"))(
      object = object,
      ...
    )
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
