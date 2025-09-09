RunCellRank <- function(
    srt = NULL,
    assay_x = "RNA",
    layer_x = "counts",
    assay_y = c("spliced", "unspliced"),
    layer_y = "counts",
    adata = NULL,
    group_by = NULL,
    n_jobs = 1,
    linear_reduction = NULL,
    nonlinear_reduction = NULL,
    basis = NULL,
    mode = "stochastic",
    fitting_by = "stochastic",
    magic_impute = FALSE,
    knn = 5,
    t = 2,
    min_shared_counts = 30,
    n_pcs = 30,
    n_neighbors = 30,
    stream_smooth = NULL,
    stream_density = 2,
    arrow_size = 5,
    arrow_length = 5,
    arrow_density = 0.5,
    calculate_velocity_genes = FALSE,
    denoise = FALSE,
    kinetics = FALSE,
    palette = "Paired",
    palcolor = NULL,
    show_plot = TRUE,
    save = FALSE,
    dpi = 300,
    dirpath = "./",
    fileprefix = "",
    return_seurat = !is.null(srt)) {
  check_python("cellrank")
  if (isTRUE(magic_impute)) {
    check_python("magic-impute")
  }
  if (all(is.null(srt), is.null(adata))) {
    log_message(
      "One of {.arg srt} or {.arg adata} must be provided",
      message_type = "error"
    )
  }
  if (is.null(group_by)) {
    log_message(
      "{.arg group_by} must be provided",
      message_type = "error"
    )
  }

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

  if (is.character(mode) && length(mode) == 1) {
    mode <- list(mode)
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
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
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
  }
  groups <- py_to_r2(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_colors(
    levels(groups) %||% unique(groups),
    palette = palette,
    palcolor = palcolor
  )

  log_message("Running {.pkg CellRank} analysis...")
  scop_analysis <- reticulate::import_from_path(
    "scop_analysis",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  adata <- do.call(scop_analysis$CellRank, args)
  log_message(
    "{.pkg CellRank} analysis completed",
    message_type = "success"
  )
  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      return(srt_append(srt_raw = srt, srt_append = srt_out))
    }
  } else {
    return(adata)
  }
}
