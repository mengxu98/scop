#' @title Run DoRothEA transcription factor activity inference
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param layer Assay layer used as the expression matrix.
#' @param species Species used to select bundled DoRothEA regulons. DoRothEA
#' only provides human and mouse regulons. For other input species, set
#' `input_species` and project expression values to this regulon species through
#' homologous gene conversion before activity inference.
#' @param input_species Species of the input expression features. If `NULL`,
#' the input is assumed to use the same gene namespace as `species`. When this
#' differs from `species`, expression features are converted with
#' [ConvertHomologs] before DoRothEA activity inference.
#' @param geneID_from_IDtype,geneID_to_IDtype Gene identifier types passed to
#' [ConvertHomologs] for cross-species projection. For bundled DoRothEA
#' regulons, `geneID_to_IDtype` should normally remain `"symbol"`.
#' @param homolog_params Additional named arguments passed to [ConvertHomologs]
#' when `input_species` differs from `species`, such as `Ensembl_version`,
#' `biomart`, `mirror`, `max_tries`, `multi_mapping`, and `collapse_fun`.
#' @param confidence DoRothEA confidence levels to keep.
#' @param regulons Optional regulon table with `tf`, `target`, `mor`, and
#' `confidence` columns. If `NULL`, bundled `dorothea_hs` or `dorothea_mm`
#' data are loaded from the `dorothea` package.
#' @param method Activity inference backend from `decoupleR`.
#' @param minsize Minimum regulon size passed to `decoupleR`.
#' @param options Additional named options passed to the selected `decoupleR`
#' function.
#' @param assay_name Name of the assay used to store TF activity scores.
#' @param new_assay Whether to store TF activity scores as a new assay.
#' @param add_meta Whether to also write TF activity scores to `srt@meta.data`
#' with the `assay_name` prefix for direct plotting with [FeatureDimPlot()].
#'
#' @return A `Seurat` object with DoRothEA results stored in
#' `srt@tools[["Dorothea"]]`, optionally TF activity scores stored in
#' `srt@meta.data`, and optionally a TF activity assay when
#' `new_assay = TRUE`. For cross-species runs, the homolog projection summary
#' is stored in `srt@tools[["Dorothea"]]$homolog_conversion`.
#' @export
#'
#' @references
#' Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D.,
#' and Saez-Rodriguez, J. (2019). Benchmark and integration of
#' resources for the estimation of human transcription factor activities.
#' \emph{Genome Research}, 29, 1363-1375. \doi{10.1101/gr.240663.118}
#'
#' Badia-i-Mompel, P., Velez Santiago, J., Braunger, J., Geiss, C.,
#' Dimitrov, D., Muller-Dott, S., Taus, P., Dugourd, A., Holland, C.H.,
#' Ramirez Flores, R.O., and Saez-Rodriguez, J. (2022). decoupleR:
#' ensemble of computational methods to infer biological activities from
#' omics data. \emph{Bioinformatics Advances}, 2, vbac016.
#' \doi{10.1093/bioadv/vbac016}
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(
#'   pancreas_sub,
#'   verbose = FALSE
#' )
#'
#' pancreas_sub <- RunDorothea(
#'   pancreas_sub,
#'   layer = "counts",
#'   species = "Mus_musculus",
#'   confidence = c("A", "B", "C"),
#'   method = "ulm",
#'   minsize = 5,
#'   new_assay = FALSE
#' )
#'
#' pancreas_sub@tools$Dorothea$regulon_summary
#' head(pancreas_sub@tools$Dorothea$result)
#'
#' activity_cols <- head(
#'   grep("^dorothea_", colnames(pancreas_sub@meta.data), value = TRUE),
#'   2
#' )
#' head(pancreas_sub@meta.data[, activity_cols, drop = FALSE])
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = activity_cols,
#'   reduction = "StandardUMAP2D",
#'   ncol = 2
#' )
RunDorothea <- function(
  srt,
  assay = NULL,
  layer = "data",
  species = c("Homo_sapiens", "Mus_musculus"),
  input_species = NULL,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  homolog_params = list(),
  confidence = c("A", "B", "C"),
  regulons = NULL,
  method = c("ulm", "viper", "wmean"),
  minsize = 5,
  options = list(),
  assay_name = "dorothea",
  new_assay = TRUE,
  add_meta = TRUE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  method <- match.arg(method)
  species <- match.arg(species)
  input_species <- input_species %||% species
  if (
    !is.character(input_species) ||
      length(input_species) != 1L ||
      is.na(input_species)
  ) {
    log_message(
      "{.arg input_species} must be a single species name",
      message_type = "error"
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)

  check_r("decoupleR", verbose = FALSE)
  if (is.null(regulons)) {
    check_r("dorothea", verbose = FALSE)
    data_name <- switch(species,
      Homo_sapiens = "dorothea_hs",
      Mus_musculus = "dorothea_mm"
    )
    env <- new.env(parent = emptyenv())
    utils::data(list = data_name, package = "dorothea", envir = env)
    regulons <- get(data_name, envir = env)
  } else {
    regulons <- as.data.frame(regulons)
  }
  if (!is.null(confidence) && "confidence" %in% colnames(regulons)) {
    regulons <- regulons[
      regulons[["confidence"]] %in% confidence, ,
      drop = FALSE
    ]
  }
  required <- c("tf", "target", "mor")
  missing <- setdiff(required, colnames(regulons))
  if (length(missing) > 0) {
    log_message(
      "{.arg regulons} must contain columns: {.val {required}}",
      message_type = "error"
    )
  }
  if (nrow(regulons) == 0) {
    log_message(
      "No DoRothEA regulon edges remain after filtering",
      message_type = "error"
    )
  }

  expr <- GetAssayData5(srt, layer = layer, assay = assay)
  expr <- expr[Matrix::rowSums(expr != 0) > 0, , drop = FALSE]
  homolog_conversion <- NULL
  if (!identical(input_species, species)) {
    invalid_homolog_params <- length(homolog_params) > 0L &&
      (is.null(names(homolog_params)) || any(!nzchar(names(homolog_params))))
    if (!is.list(homolog_params) || invalid_homolog_params) {
      log_message(
        "{.arg homolog_params} must be a named list",
        message_type = "error"
      )
    }
    reserved <- c(
      "object",
      "species_from",
      "species_to",
      "geneID_from_IDtype",
      "geneID_to_IDtype",
      "keep_unmapped",
      "verbose"
    )
    duplicated_params <- intersect(names(homolog_params), reserved)
    if (length(duplicated_params) > 0L) {
      log_message(
        "{.arg homolog_params} cannot override {.val {duplicated_params}}",
        message_type = "error"
      )
    }
    log_message(
      "Project expression features from {.val {input_species}} to {.val {species}} homologs for {.pkg DoRothEA}",
      verbose = verbose
    )
    expr <- do.call(
      ConvertHomologs.default,
      c(
        list(
          object = expr,
          species_from = input_species,
          species_to = species,
          geneID_from_IDtype = geneID_from_IDtype,
          geneID_to_IDtype = geneID_to_IDtype,
          keep_unmapped = FALSE,
          verbose = verbose
        ),
        homolog_params
      )
    )
    conversion <- attr(expr, "ConvertHomologs")
    if (!is.null(conversion)) {
      homolog_conversion <- data.frame(
        species_from = conversion$species_from,
        species_to = conversion$species_to,
        geneID_from_IDtype = paste(
          conversion$geneID_from_IDtype,
          collapse = ","
        ),
        geneID_to_IDtype = paste(
          conversion$geneID_to_IDtype,
          collapse = ","
        ),
        n_mapped_source_genes = length(unique(conversion$mapping$from_geneID)),
        n_target_homologs = length(unique(conversion$mapping$to_geneID)),
        n_unmapped_source_genes = length(conversion$unmapped),
        Ensembl_version = conversion$Ensembl_version %||% NA_character_,
        stringsAsFactors = FALSE
      )
    }
    expr <- expr[Matrix::rowSums(expr != 0) > 0, , drop = FALSE]
  }
  expr <- as.matrix(expr)
  if (nrow(expr) == 0 || ncol(expr) == 0) {
    log_message(
      "No expression values available for DoRothEA activity inference",
      message_type = "error"
    )
  }

  log_message(
    "Run {.pkg DoRothEA}/{.pkg decoupleR} with {.val {nrow(regulons)}} regulon edges",
    verbose = verbose
  )

  run_fun <- dorothea_get_run_fun(method)
  params <- c(
    list(
      mat = expr,
      network = regulons,
      .source = "tf",
      .target = "target",
      .mor = "mor",
      minsize = minsize
    ),
    options
  )
  res <- do.call(run_fun, params)
  res_df <- as.data.frame(res)
  source_col <- intersect(c("source", "tf"), colnames(res_df))[1]
  condition_col <- intersect(
    c("condition", "sample", "cell"),
    colnames(res_df)
  )[1]
  score_col <- intersect(c("score", "activity", "nes"), colnames(res_df))[1]
  if (any(is.na(c(source_col, condition_col, score_col)))) {
    log_message(
      "Unable to parse {.pkg decoupleR} result columns for DoRothEA scores",
      message_type = "error"
    )
  }
  sources <- unique(as.character(res_df[[source_col]]))
  conditions <- unique(as.character(res_df[[condition_col]]))
  scores <- matrix(
    NA_real_,
    nrow = length(sources),
    ncol = length(conditions),
    dimnames = list(sources, conditions)
  )
  idx <- cbind(
    match(as.character(res_df[[source_col]]), sources),
    match(as.character(res_df[[condition_col]]), conditions)
  )
  scores[idx] <- as.numeric(res_df[[score_col]])
  missing_cells <- setdiff(colnames(srt), colnames(scores))
  if (length(missing_cells) > 0L) {
    log_message(
      "{.pkg decoupleR} did not return scores for all cells in {.arg srt}",
      message_type = "error"
    )
  }
  scores <- scores[, colnames(srt), drop = FALSE]

  if (isTRUE(new_assay)) {
    srt[[assay_name]] <- Seurat::CreateAssayObject(data = scores)
    srt[[assay_name]] <- Seurat::AddMetaData(
      object = srt[[assay_name]],
      metadata = data.frame(
        termnames = rownames(scores),
        row.names = rownames(scores),
        stringsAsFactors = FALSE
      )
    )
    log_message(
      "{.pkg DoRothEA} TF activity scores stored in assay {.val {assay_name}}",
      verbose = verbose
    )
  }
  if (isTRUE(add_meta)) {
    meta_scores <- as.data.frame(t(scores), check.names = FALSE)
    colnames(meta_scores) <- make.names(
      paste(assay_name, colnames(meta_scores), sep = "_")
    )
    srt <- Seurat::AddMetaData(srt, metadata = meta_scores)
    log_message(
      "{.pkg DoRothEA} TF activity scores stored in {.cls Seurat} metadata",
      verbose = verbose
    )
  }

  srt@tools[["Dorothea"]] <- list(
    scores = scores,
    result = res_df,
    regulon_summary = data.frame(
      n_tfs = length(unique(regulons[["tf"]])),
      n_targets = length(unique(regulons[["target"]])),
      n_edges = nrow(regulons),
      confidence = if ("confidence" %in% colnames(regulons)) {
        paste(sort(unique(regulons[["confidence"]])), collapse = ",")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    ),
    parameters = list(
      assay = assay,
      layer = layer,
      species = species,
      input_species = input_species,
      geneID_from_IDtype = geneID_from_IDtype,
      geneID_to_IDtype = geneID_to_IDtype,
      confidence = confidence,
      method = method,
      minsize = minsize,
      assay_name = assay_name,
      new_assay = new_assay,
      add_meta = add_meta,
      homolog_params = homolog_params,
      options = options
    ),
    homolog_conversion = homolog_conversion
  )
  srt
}

dorothea_get_run_fun <- function(method) {
  switch(method,
    ulm = getExportedValue("decoupleR", "run_ulm"),
    viper = getExportedValue("decoupleR", "run_viper"),
    wmean = getExportedValue("decoupleR", "run_wmean")
  )
}

#' @title Plot differential DoRothEA TF activity
#'
#' @description
#' Compare DoRothEA transcription factor activity between two groups and draw
#' a signed bar plot. Bar height is the mean activity difference
#' `group1 - group2`; fill color is the signed `-log10(p)` or
#' `-log10(adjusted p)`, where the sign follows the activity difference.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object containing results from [RunDorothea()].
#' @param group.by Metadata column used to define groups.
#' @param group1,group2 Two group labels to compare. Positive logFC means
#' higher TF activity in `group1`.
#' @param tool_name Name of the `srt@tools` entry created by [RunDorothea()].
#' @param features TFs to plot. If `NULL`, all TFs are tested and the top
#' `top_n` TFs are shown.
#' @param top_n Number of TFs to show when `features = NULL`. Set `NULL` to
#' show all tested TFs.
#' @param test.use Statistical test used for each TF.
#' @param p.adjust.method Method passed to [stats::p.adjust].
#' @param color.by P-value column used for bar fill.
#' @param rank.by Metric used to select top TFs when `features = NULL`.
#' @param sort.by Metric used to order TFs on the x-axis.
#' @param p_floor Lower bound used before `-log10()` transformation.
#' @param bar_width Width of bars.
#' @param cols Color vector of length 3 for low, midpoint, and high values.
#' @param title,xlab,ylab,fill.title Axis, plot, and legend titles.
#' @param angle X-axis text angle.
#' @param hjust,vjust X-axis text justification.
#' @param theme_use Theme function used to style the plot.
#' @param theme_args Other arguments passed to `theme_use`.
#' @param return_data Whether to return a list with the plot and statistics.
#'
#' @return A `ggplot` object, or a list with `plot` and `data` when
#' `return_data = TRUE`.
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub, verbose = FALSE)
#' pancreas_sub <- RunDorothea(
#'   pancreas_sub,
#'   layer = "counts",
#'   species = "Mus_musculus",
#'   method = "ulm",
#'   minsize = 5,
#'   new_assay = FALSE
#' )
#' groups <- unique(as.character(pancreas_sub$CellType))
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group1 = groups[1],
#'   group2 = groups[2]
#' )
DorotheaPlot <- function(
  srt,
  group.by,
  group1,
  group2,
  tool_name = "Dorothea",
  features = NULL,
  top_n = 30,
  test.use = c("wilcox.test", "t.test"),
  p.adjust.method = "BH",
  color.by = c("p_val", "p_val_adj"),
  rank.by = c("abs_logFC", "p_val", "p_val_adj", "logFC"),
  sort.by = c("logFC", "abs_logFC", "p_val", "p_val_adj"),
  p_floor = .Machine$double.xmin,
  bar_width = 0.85,
  cols = c("#2166AC", "white", "#B2182B"),
  title = NULL,
  xlab = NULL,
  ylab = "logFC",
  fill.title = NULL,
  angle = 90,
  hjust = 1,
  vjust = 0.5,
  theme_use = "theme_scop",
  theme_args = list(),
  return_data = FALSE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (length(group.by) != 1L || !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} must be one metadata column in {.arg srt}",
      message_type = "error"
    )
  }
  if (
    !tool_name %in% names(srt@tools) || is.null(srt@tools[[tool_name]]$scores)
  ) {
    log_message(
      "No DoRothEA scores found in {.code srt@tools[[{tool_name}]]$scores}",
      message_type = "error"
    )
  }
  test.use <- match.arg(test.use)
  color.by <- match.arg(color.by)
  rank.by <- match.arg(rank.by)
  sort.by <- match.arg(sort.by)
  if (
    length(group1) != 1L ||
      length(group2) != 1L ||
      any(is.na(c(group1, group2)))
  ) {
    log_message(
      "{.arg group1} and {.arg group2} must each be a single group label",
      message_type = "error"
    )
  }
  if (!is.null(top_n)) {
    if (length(top_n) != 1L || is.na(top_n) || top_n < 1) {
      log_message(
        "{.arg top_n} must be a positive integer or {.code NULL}",
        message_type = "error"
      )
    }
    top_n <- as.integer(top_n)
  }
  if (length(cols) != 3L) {
    log_message(
      "{.arg cols} must contain three colors: low, midpoint, and high",
      message_type = "error"
    )
  }

  scores <- as.matrix(srt@tools[[tool_name]]$scores)
  groups <- as.character(srt@meta.data[[group.by]])
  names(groups) <- rownames(srt@meta.data)
  cells <- intersect(colnames(scores), names(groups))
  cells <- cells[groups[cells] %in% c(group1, group2)]
  if (length(cells) == 0L) {
    log_message(
      "No cells from {.val {group1}} or {.val {group2}} are shared between DoRothEA scores and metadata",
      message_type = "error"
    )
  }
  cells1 <- cells[groups[cells] == group1]
  cells2 <- cells[groups[cells] == group2]
  if (length(cells1) == 0L || length(cells2) == 0L) {
    log_message(
      "Both {.arg group1} and {.arg group2} must contain cells",
      message_type = "error"
    )
  }

  features_is_null <- is.null(features)
  if (features_is_null) {
    features <- rownames(scores)
  } else {
    features <- unique(as.character(features))
    missing_features <- setdiff(features, rownames(scores))
    if (length(missing_features) > 0L) {
      log_message(
        "Dropping TFs not found in DoRothEA scores: {.val {missing_features}}",
        message_type = "warning",
        verbose = verbose
      )
    }
    features <- intersect(features, rownames(scores))
  }
  if (length(features) == 0L) {
    log_message(
      "No TFs are available for plotting",
      message_type = "error"
    )
  }

  log_message(
    "Compare DoRothEA TF activity: {.val {group1}} vs {.val {group2}}",
    verbose = verbose
  )
  mat1 <- scores[features, cells1, drop = FALSE]
  mat2 <- scores[features, cells2, drop = FALSE]
  mean1 <- Matrix::rowMeans(mat1, na.rm = TRUE)
  mean2 <- Matrix::rowMeans(mat2, na.rm = TRUE)
  p_val <- vapply(
    features,
    function(tf) {
      x <- as.numeric(mat1[tf, ])
      y <- as.numeric(mat2[tf, ])
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
      if (length(x) < 1L || length(y) < 1L) {
        return(NA_real_)
      }
      tryCatch(
        switch(test.use,
          wilcox.test = stats::wilcox.test(x, y)$p.value,
          t.test = stats::t.test(x, y)$p.value
        ),
        error = function(e) NA_real_
      )
    },
    numeric(1)
  )
  stat_df <- data.frame(
    TF = features,
    group1 = group1,
    group2 = group2,
    mean1 = as.numeric(mean1[features]),
    mean2 = as.numeric(mean2[features]),
    logFC = as.numeric(mean1[features] - mean2[features]),
    p_val = p_val,
    stringsAsFactors = FALSE
  )
  stat_df$p_val[!is.finite(stat_df$p_val) | stat_df$p_val < 0] <- NA_real_
  stat_df$p_val_adj <- stats::p.adjust(stat_df$p_val, method = p.adjust.method)
  stat_df$neglog10_p_val <- -log10(pmax(stat_df$p_val, p_floor, na.rm = TRUE))
  stat_df$neglog10_p_val_adj <- -log10(pmax(
    stat_df$p_val_adj,
    p_floor,
    na.rm = TRUE
  ))
  stat_df$signed_neglog10_p_val <- sign(stat_df$logFC) * stat_df$neglog10_p_val
  stat_df$signed_neglog10_p_val_adj <- sign(stat_df$logFC) *
    stat_df$neglog10_p_val_adj
  stat_df$neglog10_p_val[!is.finite(stat_df$neglog10_p_val)] <- 0
  stat_df$neglog10_p_val_adj[!is.finite(stat_df$neglog10_p_val_adj)] <- 0
  stat_df$signed_neglog10_p_val[!is.finite(stat_df$signed_neglog10_p_val)] <- 0
  stat_df$signed_neglog10_p_val_adj[
    !is.finite(stat_df$signed_neglog10_p_val_adj)
  ] <- 0

  if (isTRUE(features_is_null) && !is.null(top_n)) {
    top_n <- min(as.integer(top_n), nrow(stat_df))
    if (rank.by == "abs_logFC") {
      rank_order <- order(-abs(stat_df$logFC), stat_df$p_val, na.last = TRUE)
    } else if (rank.by == "logFC") {
      rank_order <- order(-stat_df$logFC, stat_df$p_val, na.last = TRUE)
    } else {
      rank_order <- order(
        stat_df[[rank.by]],
        -abs(stat_df$logFC),
        na.last = TRUE
      )
    }
    stat_df <- stat_df[rank_order[seq_len(top_n)], , drop = FALSE]
  }

  if (sort.by == "abs_logFC") {
    plot_order <- order(-abs(stat_df$logFC), -stat_df$logFC, na.last = TRUE)
  } else if (sort.by == "logFC") {
    plot_order <- order(stat_df$logFC, decreasing = TRUE)
  } else {
    plot_order <- order(stat_df[[sort.by]], -abs(stat_df$logFC), na.last = TRUE)
  }
  stat_df <- stat_df[plot_order, , drop = FALSE]
  stat_df$TF <- factor(stat_df$TF, levels = stat_df$TF)

  fill_col <- switch(color.by,
    p_val = "signed_neglog10_p_val",
    p_val_adj = "signed_neglog10_p_val_adj"
  )
  max_fill <- max(abs(stat_df[[fill_col]]), na.rm = TRUE)
  if (!is.finite(max_fill) || max_fill == 0) {
    max_fill <- 1
  }
  title <- title %||% paste(group1, "vs.", group2)
  fill.title <- fill.title %||%
    ifelse(
      color.by == "p_val",
      "-log10(p)",
      "-log10(padj)"
    )
  p <- ggplot2::ggplot(
    stat_df,
    ggplot2::aes(x = TF, y = logFC, fill = .data[[fill_col]])
  ) +
    ggplot2::geom_col(width = bar_width) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(
      low = cols[1],
      mid = cols[2],
      high = cols[3],
      midpoint = 0,
      limits = c(-max_fill, max_fill),
      labels = abs,
      name = fill.title
    ) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    do.call(theme_use, theme_args) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = angle,
        hjust = hjust,
        vjust = vjust
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )

  if (isTRUE(return_data)) {
    return(list(plot = p, data = stat_df))
  }
  p
}
