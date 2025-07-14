#' @title Run Palantir analysis
#'
#' @inheritParams RunPAGA
#' @param dm_n_components The number of diffusion components to calculate.
#' @param dm_alpha Normalization parameter for the diffusion operator.
#' @param dm_n_eigs Number of eigen vectors to use.
#' @param early_group Name of the group to start Palantir analysis from.
#' @param early_cell Name of the cell to start Palantir analysis from.
#' @param terminal_groups Character vector specifying terminal groups for Palantir analysis.
#' @param terminal_cells Character vector specifying terminal cells for Palantir analysis.
#' @param num_waypoints Number of waypoints to be included.
#' @param scale_components Should the cell fate probabilities be scaled for each component independently?
#' @param use_early_cell_as_start Should the starting cell for each terminal group be set as early_cell?
#' @param adjust_early_cell Whether to adjust the early cell to the cell with the minimum pseudotime value.
#' @param adjust_terminal_cells Whether to adjust the terminal cells to the cells with the maximum pseudotime value for each terminal group.
#' @param max_iterations Maximum number of iterations for pseudotime convergence.
#' @param n_jobs The number of parallel jobs to run.
#'
#' @seealso \code{\link{srt_to_adata}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunPalantir(
#'   srt = pancreas_sub,
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   early_group = "Ductal",
#'   use_early_cell_as_start = TRUE,
#'   terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon")
#' )
#' head(pancreas_sub[[]])
#' FeatureDimPlot(
#'   pancreas_sub,
#'   c("palantir_pseudotime", "palantir_diff_potential")
#' )
#' FeatureDimPlot(
#'   pancreas_sub,
#'   paste0(
#'     c("Alpha", "Beta", "Delta", "Epsilon"),
#'     "_diff_potential"
#'   )
#' )
#' }
RunPalantir <- function(
    srt = NULL,
    assay_x = "RNA",
    layer_x = "counts",
    assay_y = c("spliced", "unspliced"),
    layer_y = "counts",
    adata = NULL,
    group_by = NULL,
    linear_reduction = NULL,
    nonlinear_reduction = NULL,
    basis = NULL,
    n_pcs = 30,
    n_neighbors = 30,
    dm_n_components = 10,
    dm_alpha = 0,
    dm_n_eigs = NULL,
    early_group = NULL,
    early_cell = NULL,
    terminal_cells = NULL,
    terminal_groups = NULL,
    num_waypoints = 1200,
    scale_components = TRUE,
    use_early_cell_as_start = TRUE,
    adjust_early_cell = FALSE,
    adjust_terminal_cells = FALSE,
    max_iterations = 25,
    n_jobs = 1,
    point_size = 20,
    palette = "Paired",
    palcolor = NULL,
    show_plot = TRUE,
    save = FALSE,
    dpi = 300,
    dirpath = "./",
    fileprefix = "",
    return_seurat = !is.null(srt)) {
  check_python("palantir")
  if (all(is.null(srt), is.null(adata))) {
    log_message(
      "{.arg srt} or {.arg adata} must be provided.",
      message_type = "error"
    )
  }
  if (is.null(group_by) && any(!is.null(early_group), !is.null(terminal_groups))) {
    log_message(
      "{.arg group_by} must be provided when {.arg early_group} or {.arg terminal_groups} provided.",
      message_type = "error"
    )
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    log_message(
      "{.arg linear_reduction} or {.arg nonlinear_reduction} must be provided at least one.",
      message_type = "error"
    )
  }
  if (is.null(early_cell) && is.null(early_group)) {
    log_message(
      "{.arg early_cell} or {.arg early_group} must be provided.",
      message_type = "error"
    )
  }
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call_envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call_envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call_envir)
    } else {
      arg
    }
  })
  args <- args[
    !names(args) %in%
      c(
        "srt",
        "assay_x",
        "layer_x",
        "assay_y",
        "layer_y",
        "return_seurat",
        "palette",
        "palcolor"
      )
  ]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_x = assay_x,
      layer_x = layer_x,
      assay_y = assay_y,
      layer_y = layer_y
    )
    if (!is.null(linear_reduction)) {
      if (!startsWith(linear_reduction, "X_")) {
        args[["linear_reduction"]] <- paste0(
          "X_", tolower(linear_reduction)
        )
      }
    }
    if (!is.null(nonlinear_reduction)) {
      if (!startsWith(nonlinear_reduction, "X_")) {
        args[["nonlinear_reduction"]] <- paste0(
          "X_", tolower(nonlinear_reduction)
        )
      }
    }
    if (is.null(basis)) {
      if (!is.null(args[["nonlinear_reduction"]])) {
        args[["basis"]] <- args[["nonlinear_reduction"]]
      } else if (!is.null(args[["linear_reduction"]])) {
        args[["basis"]] <- args[["linear_reduction"]]
      }
    } else {
      if (!startsWith(basis, "X_")) {
        args[["basis"]] <- paste0("X_", tolower(basis))
      }
    }
  }
  groups <- py_to_r2(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scop(
    levels(groups) %||% unique(groups),
    palette = palette,
    palcolor = palcolor
  )

  scop_analysis <- reticulate::import_from_path(
    "scop_analysis",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  adata <- do.call(scop_analysis$Palantir, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- srt_append(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- srt_append(
        srt_raw = srt_out1,
        srt_append = srt_out,
        pattern = "(palantir)|(dm_kernel)|(_diff_potential)",
        overwrite = TRUE,
        verbose = FALSE
      )
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}
