#' @title Run Palantir analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunCellRank
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
#' @param point_size The point size for plotting.
#'
#' @seealso [srt_to_adata]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunPalantir(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   early_group = "Ductal",
#'   terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon")
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   c("palantir_pseudotime", "palantir_diff_potential")
#' )
#'
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
    group.by = NULL,
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
    cores = 1,
    point_size = 20,
    palette = "Paired",
    palcolor = NULL,
    legend.position = "on data",
    show_plot = FALSE,
    save_plot = FALSE,
    plot_format = c("pdf", "png", "svg"),
    plot_dpi = 300,
    plot_prefix = "palantir",
    dirpath = "./",
    return_seurat = !is.null(srt),
    verbose = TRUE) {
  PrepareEnv()
  check_python("palantir", verbose = verbose)
  if (all(is.null(srt), is.null(adata))) {
    log_message(
      "{.arg srt} or {.arg adata} must be provided.",
      message_type = "error"
    )
  }
  if (is.null(group.by) && any(!is.null(early_group), !is.null(terminal_groups))) {
    log_message(
      "{.arg group.by} must be provided when {.arg early_group} or {.arg terminal_groups} provided.",
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

  args[["legend_loc"]] <- legend.position

  args[["n_jobs"]] <- as.integer(cores)

  args[["save"]] <- save_plot
  args[["dpi"]] <- plot_dpi
  args[["fileprefix"]] <- plot_prefix
  params <- c(
    "srt",
    "assay_x",
    "layer_x",
    "assay_y",
    "layer_y",
    "return_seurat",
    "palette",
    "palcolor",
    "save_plot",
    "plot_format",
    "plot_dpi",
    "plot_prefix",
    "legend.position",
    "cores"
  )
  args <- args[!names(args) %in% params]

  if (!is.null(srt)) {
    if (is.null(linear_reduction)) {
      linear_reduction <- DefaultReduction(srt)
    } else {
      linear_reduction <- DefaultReduction(srt, pattern = linear_reduction)
    }
    if (!linear_reduction %in% names(srt@reductions)) {
      log_message(
        "{.val {linear_reduction}} is not in the srt reduction names",
        message_type = "error"
      )
    }

    if (is.null(nonlinear_reduction)) {
      nonlinear_reduction <- DefaultReduction(srt)
    } else {
      nonlinear_reduction <- DefaultReduction(srt, pattern = nonlinear_reduction)
    }
    if (!nonlinear_reduction %in% names(srt@reductions)) {
      log_message(
        "{.val {nonlinear_reduction}} is not in the srt reduction names",
        message_type = "error"
      )
    }

    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_x = assay_x,
      layer_x = layer_x,
      assay_y = assay_y,
      layer_y = layer_y
    )

    if (!is.null(linear_reduction)) {
      args[["linear_reduction"]] <- linear_reduction
    }
    if (!is.null(nonlinear_reduction)) {
      args[["nonlinear_reduction"]] <- nonlinear_reduction
    }
    if (is.null(basis)) {
      if (!is.null(nonlinear_reduction)) {
        args[["basis"]] <- nonlinear_reduction
      } else if (!is.null(linear_reduction)) {
        args[["basis"]] <- linear_reduction
      }
    } else {
      args[["basis"]] <- basis
    }
  }
  if ("group.by" %in% names(args)) {
    args[["group_by"]] <- args[["group.by"]]
    args[["group.by"]] <- NULL
  }
  groups <- py_to_r2(args[["adata"]]$obs)[[group.by]]
  args[["palette"]] <- palette_colors(
    levels(groups) %||% unique(groups),
    palette = palette,
    palcolor = palcolor
  )

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  adata <- do.call(functions$Palantir, args)

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
