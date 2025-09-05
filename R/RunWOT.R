#' Run WOT analysis
#'
#' @inheritParams RunPAGA
#' @param time_field A character string specifying the column name in `adata.obs` or `srt@meta.data` that contains the time information.
#' @param growth_iters An integer specifying the number of growth iterations to perform during the OT Model computation. Default is 3.
#' @param tmap_out A character string specifying the path to store the computed transport maps.
#' @param time_from A numeric value specifying the starting time point for trajectory and fate analysis.
#' @param time_to A numeric value specifying the ending time point for trajectory and fate analysis. If not provided, only trajectory and fate analysis for the specified `time_from` will be performed.
#' @param get_coupling A logical value indicating whether to compute and store the coupling matrix between the specified `time_from` and `time_to`. Default is FALSE.
#' @param recalculate A logical value indicating whether to recalculate the transport maps even if they already exist at the specified `tmap_out` location. Default is FALSE.
#'
#' @seealso \link{srt_to_adata}
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#'
#' print(range(pancreas_sub$Lineage1, na.rm = TRUE))
#'
#' pancreas_sub <- RunWOT(
#'   pancreas_sub,
#'   group_by = "SubCellType",
#'   time_field = "Lineage1",
#'   time_from = min(pancreas_sub$Lineage1, na.rm = TRUE),
#'   time_to = max(pancreas_sub$Lineage1, na.rm = TRUE),
#'   get_coupling = TRUE,
#'   tmap_out = "tmaps/lineage_tmap"
#' )
#'
#' pancreas_sub$Custom_Time <- sample(
#'   1:10,
#'   ncol(pancreas_sub),
#'   replace = TRUE
#' )
#' pancreas_sub <- RunWOT(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   time_field = "Custom_Time",
#'   time_from = 1,
#'   time_to = 10,
#'   tmap_out = "tmaps/custom_tmap"
#' )
#' }
RunWOT <- function(
    srt = NULL,
    assay_x = "RNA",
    layer_x = "counts",
    assay_y = c("spliced", "unspliced"),
    layer_y = "counts",
    adata = NULL,
    group_by = NULL,
    time_field = "Time",
    growth_iters = 3L,
    tmap_out = "tmaps/tmap_out",
    time_from = NULL,
    time_to = NULL,
    get_coupling = FALSE,
    recalculate = FALSE,
    palette = "Paired",
    palcolor = NULL,
    show_plot = TRUE,
    save = FALSE,
    dpi = 300,
    dirpath = "./",
    fileprefix = "",
    return_seurat = !is.null(srt)) {
  check_python("wot")
  if (all(is.null(srt), is.null(adata))) {
    log_message(
      "One of 'srt', 'adata' must be provided.",
      message_type = "error"
    )
  }
  if (is.null(group_by)) {
    log_message(
      "'group_by' must be provided.",
      message_type = "error"
    )
  }
  if (is.null(time_field)) {
    log_message(
      "'time_field' must be provided.",
      message_type = "error"
    )
  }
  if (is.null(time_from)) {
    log_message(
      "'time_from' must be provided.",
      message_type = "error"
    )
  }
  if (isTRUE(get_coupling) && is.null(time_to)) {
    log_message(
      "The 'get_coupling' paramter is only valid when 'time_to' is specified.",
      message_type = "warning"
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

  scop_analysis <- reticulate::import_from_path(
    "scop_analysis",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  adata <- do.call(scop_analysis$WOT, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- srt_append(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- srt_append(
        srt_raw = srt_out1,
        srt_append = srt_out,
        pattern = "(trajectory_)|(fates_)|(transition_)|(coupling_)",
        overwrite = TRUE,
        verbose = FALSE
      )
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}
