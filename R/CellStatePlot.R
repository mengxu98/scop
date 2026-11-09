#' @title Composite cell-state plot
#'
#' @description
#' Combine a cluster similarity tree, class and subclass annotations, program
#' assignment, sample contributions, marker expression and precomputed GSEA.
#' Columns represent clusters throughout the panel. The layout uses compact
#' physical spacing and preserves its aspect ratio when resized.
#'
#' @md
#' @param object An analysed Seurat or Python AnnData object. Read saved objects
#'   first with `readRDS()` or `h5ad_to_srt()`. AnnData uses the existing Python
#'   environment; this function installs nothing.
#' @param group.by,sample.by Metadata columns identifying clusters and biological
#'   samples or donors. Both must be complete.
#' @param features Marker genes in display order. Unmatched genes remain missing.
#' @param program.by A metadata column containing per-cell program assignments.
#' @param program.scores Numeric metadata columns containing program scores.
#'   Names, if supplied, are used as program labels. The highest score wins;
#'   ties are Ambiguous and cells with no finite scores are Unassigned.
#' @param programs Alternatively, a uniquely named list of gene sets. Each gene
#'   is standardized across cells using the population SD, then averaged within
#'   each set. The largest mean determines the relative assignment. At least
#'   two sets must have min.program.genes matched, variable genes.
#'   Supply exactly one of program.by, program.scores and programs.
#' @param gsea A precomputed RunGSEA result or enrichment table. Native columns
#'   are Groups, Database, ID, Description, NES, pvalue and p.adjust.
#'   Tidy columns cluster, feature, nes, p_value and optional padj, status
#'   are also accepted. NOT_EVALUABLE records remain missing.
#'   If NULL, use the Seurat object's stored RunGSEA result for group.by and
#'   test.use. This plot does not run differential expression or enrichment.
#' @param pathways Pathway IDs or unambiguous descriptions in display order.
#'   NULL displays all supplied pathways, without a significance filter.
#' @param db Database to select when the GSEA table contains several databases.
#' @param test.use Test used to locate stored RunGSEA results.
#' @param class.by,subclass.by Optional metadata columns giving a single
#'   annotation per cluster. Mixed annotations require an explicit map.
#' @param class.map,subclass.map Optional named character vectors mapping
#'   cluster IDs to annotations. These take precedence over metadata columns.
#' @param assay Seurat assay. NULL uses DefaultAssay(object).
#' @param layer Seurat expression layer (default data), or AnnData layer
#'   (default X). Select a single joined layer explicitly.
#' @param input.scale Either log1p normalized, nonnegative expression, or counts.
#'   Counts are normalized to 10,000 per cell and transformed with log1p.
#'   Scaled expression, residuals and integrated values are inappropriate.
#' @param reduction Existing Seurat reduction or AnnData obsm key for the tree.
#'   NULL uses up to 500 variable cluster-average gene profiles.
#' @param tree Optional hclust with unique cluster-ID labels. Otherwise use
#'   Pearson distance and complete linkage on cluster-average profiles.
#'   This is a similarity tree, not evidence of time or lineage.
#' @param min.pct Minimum expressing fraction for a marker dot.
#' @param expression.cutoff Expression strictly above this value counts as
#'   expressed, after any counts normalization.
#' @param min.program.genes Minimum matched variable genes per program.
#' @param feature.labels,pathway.labels Optional named vectors mapping feature
#'   or pathway IDs to display labels. Newlines are supported for pathways.
#' @param palette Categorical palette name passed to `thisplot::palette_colors()`.
#' @param expression_palette,enrichment_palette Continuous palette names for
#'   expression and NES. The default RdBu shows negative NES in blue and
#'   positive NES in red.
#' @param palcolor Optional list of named color vectors for class, subclass,
#'   program and sample, or color ramps for expression and enrichment.
#' @param family Font family. Systemfonts, when available, supplies accurate
#'   text measurements; otherwise approximate character widths are used.
#' @param title,panel_label Optional title and panel letter. Both default to NULL.
#' @param return_data Return the plot, numerical summaries, tree, provenance
#'   and native dimensions instead of only the plot.
#'
#' @return A printable grid grob accepted by ggplot2::ggsave(). Attributes
#'   width_mm and height_mm give its native dimensions. With return_data = TRUE,
#'   a list with plot, data, tree, provenance, width_mm and height_mm.
#'
#' @details
#' Marker colors encode per-gene z scores of cluster means, using population SD;
#' zero-variance genes have z = 0. Colors clip to 0--2. Marker area increases
#' with expressing fraction over the displayed range. GSEA colors encode NES
#' clipped to -2.5--2.5; area is proportional to -ln(raw P). P = 1 has zero
#' area. P = 0 uses 1e-10 for drawing only. Raw values are retained.
#' Missing expression or enrichment combinations use gray crosses.
#'
#' Sample and program bars are fractions of cells within a cluster; they are
#' descriptive compositions, not sample-level statistical tests. Gene-set
#' program assignments are relative scores, not inferred biological states.
#' Only the selected expression layer is read from AnnData, into R memory.
#' The input object is not modified and constructing a plot writes no files.
#'
#' @seealso [adata_to_srt()], [h5ad_to_srt()], [RunGSEA()]
#' @export
#'
#' @examples
#' counts <- rbind(
#'   M1 = c(2, 0, 0, 0, 4, 0),
#'   M2 = c(0, 0, 2, 2, 0, 0),
#'   M3 = c(0, 1, 0, 2, 1, 3),
#'   M4 = rep(1, 6)
#' )
#' colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
#' object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
#' object$cluster <- rep(c("A", "B", "C"), each = 2)
#' object$donor <- c("d1", "d2", "d1", "d2", "d1", "d1")
#' object$program <- c("P1", "P1", "P2", "P2", "P1", "P2")
#' enrichment <- expand.grid(
#'   Groups = c("A", "B", "C"), ID = c("Pathway1", "Pathway2", "Pathway3")
#' )
#' enrichment$NES <- c(1.5, -1.2, 0.5, 1, 0.7, -1.8, -0.6, 1.4, 2)
#' enrichment$pvalue <- c(0.05, 0.01, 0.2, 0.02, 0.3, 0.01, 0.1, 0.05, NA)
#' CellStatePlot(
#'   object,
#'   group.by = "cluster", sample.by = "donor",
#'   features = rownames(counts), program.by = "program",
#'   gsea = enrichment, layer = "counts", input.scale = "counts"
#' )
CellStatePlot <- function(
  object,
  group.by,
  sample.by,
  features,
  program.by = NULL,
  program.scores = NULL,
  programs = NULL,
  gsea = NULL,
  pathways = NULL,
  db = NULL,
  test.use = "wilcox",
  class.by = NULL,
  subclass.by = NULL,
  class.map = NULL,
  subclass.map = NULL,
  assay = NULL,
  layer = NULL,
  input.scale = c("log1p", "counts"),
  reduction = NULL,
  tree = NULL,
  min.pct = 0.25,
  expression.cutoff = 0,
  min.program.genes = 3,
  feature.labels = NULL,
  pathway.labels = NULL,
  palette = "Chinese",
  expression_palette = "Reds",
  enrichment_palette = "RdBu",
  palcolor = NULL,
  family = "sans",
  title = NULL,
  panel_label = NULL,
  return_data = FALSE
) {
  input.scale <- match.arg(input.scale)
  for (nm in c("group.by", "sample.by")) {
    value <- get(nm)
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
      stop(nm, " must name one metadata column.")
    }
  }
  if (!is.numeric(min.pct) || length(min.pct) != 1L ||
    !is.finite(min.pct) || min.pct < 0 || min.pct > 1) {
    stop("min.pct must be between 0 and 1.")
  }
  if (length(expression.cutoff) != 1L || !is.finite(expression.cutoff) ||
    expression.cutoff < 0) {
    stop("expression.cutoff must be nonnegative.")
  }
  if (length(min.program.genes) != 1L || !is.finite(min.program.genes) ||
    min.program.genes < 1 || min.program.genes %% 1 != 0) {
    stop("min.program.genes must be a positive integer.")
  }
  choices <- list(program.by, program.scores, programs)
  if (sum(!vapply(choices, is.null, logical(1))) != 1L) {
    stop("Supply exactly one of program.by, program.scores and programs.")
  }
  if (identical(layer, "scale.data")) {
    stop("Use a nonnegative log1p expression layer or counts, not scale.data.")
  }
  input <- cellstate_input(object, assay, layer, reduction)
  if (is.null(gsea) && !is.null(input$gsea.object)) {
    key <- paste("GSEA", group.by, test.use, sep = "_")
    gsea <- input$gsea.object@tools[[key]]
  }
  if (is.null(gsea)) {
    stop("Provide precomputed gsea results, or run RunGSEA() before plotting.")
  }
  x <- input$x
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
    anyNA(rownames(x)) || anyNA(colnames(x)) ||
    any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
    anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop("Expression needs unique, nonempty gene and cell IDs.")
  }
  meta <- input$meta
  if (!setequal(colnames(x), rownames(meta))) {
    stop("Expression and metadata cell IDs do not match.")
  }
  meta <- meta[colnames(x), , drop = FALSE]
  required <- c(
    group.by, sample.by,
    if (is.null(class.map)) class.by,
    if (is.null(subclass.map)) subclass.by,
    program.by, program.scores
  )
  if (!all(required %in% names(meta))) {
    stop("Missing metadata columns: ", paste(setdiff(required, names(meta)), collapse = ", "))
  }
  for (nm in c(group.by, sample.by)) {
    if (anyNA(meta[[nm]]) || any(!nzchar(as.character(meta[[nm]])))) {
      stop(nm, " contains missing or empty values.")
    }
  }
  x <- methods::as(methods::as(
    Matrix::Matrix(x, sparse = TRUE),
    "generalMatrix"
  ), "CsparseMatrix")
  if (any(!is.finite(x@x)) || any(x@x < 0)) {
    stop("Expression must be finite and nonnegative; do not use scaled residuals.")
  }
  if (input.scale == "counts") {
    totals <- Matrix::colSums(x)
    if (any(totals <= 0)) stop("Counts contain cells with zero library size.")
    ids <- dimnames(x)
    x <- x %*% Matrix::Diagonal(x = 1e4 / totals)
    dimnames(x) <- ids
    x@x <- log1p(x@x)
  }
  groups <- sort(unique(as.character(meta[[group.by]])))
  if (length(groups) < 2L) stop("At least two clusters are required.")
  code <- match(as.character(meta[[group.by]]), groups)
  sizes <- tabulate(code, length(groups))
  membership <- Matrix::sparseMatrix(
    i = seq_len(ncol(x)), j = code, x = 1,
    dims = c(ncol(x), length(groups))
  )
  means <- as.matrix(x %*% membership) / rep(sizes, each = nrow(x))
  dimnames(means) <- list(rownames(x), groups)
  clusters <- data.frame(
    cluster = groups,
    class = cellstate_annotation(meta, group.by, groups, class.by, class.map, "Cells"),
    subclass = cellstate_annotation(meta, group.by, groups, subclass.by, subclass.map, groups),
    n_cells = sizes,
    n_samples = vapply(groups, function(g) {
      length(unique(meta[[sample.by]][as.character(meta[[group.by]]) == g]))
    }, integer(1)),
    stringsAsFactors = FALSE
  )
  features <- unique(as.character(unlist(features, use.names = FALSE)))
  if (!length(features) || anyNA(features) || any(!nzchar(features))) {
    stop("features must contain nonempty gene IDs.")
  }
  matched <- intersect(features, rownames(x))
  if (!length(matched)) stop("No marker genes matched the expression gene IDs.")
  absent <- setdiff(features, matched)
  if (length(absent)) {
    warning("Unmatched markers retained as missing: ", paste(absent, collapse = ", "))
  }
  gm <- means[matched, , drop = FALSE]
  centered <- gm - rowMeans(gm)
  sd <- sqrt(rowMeans(centered^2))
  z <- centered / ifelse(sd > 0, sd, 1)
  fraction <- as.matrix(
    (x[matched, , drop = FALSE] > expression.cutoff) %*% membership
  ) / rep(sizes, each = length(matched))
  markers <- expand.grid(cluster = groups, feature = features, stringsAsFactors = FALSE)
  index <- match(markers$feature, matched) +
    (match(markers$cluster, groups) - 1L) * length(matched)
  markers$mean <- gm[index]
  markers$mean_z <- z[index]
  markers$proportion <- fraction[index]
  coverage <- NULL
  if (!is.null(program.by)) {
    if (length(program.by) != 1L) stop("program.by must name one metadata column.")
    assignment <- as.character(meta[[program.by]])
  } else if (!is.null(program.scores)) {
    if (!length(program.scores) ||
      !all(vapply(meta[program.scores], is.numeric, logical(1)))) {
      stop("program.scores must name numeric metadata columns.")
    }
    labels <- names(program.scores)
    if (is.null(labels)) labels <- program.scores
    if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
      stop("Program labels must be unique and nonempty.")
    }
    assignment <- cellstate_winners(as.matrix(meta[program.scores]), labels)
  } else {
    if (!is.list(programs) || is.null(names(programs)) ||
      anyNA(names(programs)) || any(!nzchar(names(programs))) ||
      anyDuplicated(names(programs))) {
      stop("programs must be a uniquely named list of gene sets.")
    }
    union <- intersect(unique(as.character(unlist(programs))), rownames(x))
    px <- x[union, , drop = FALSE]
    mu <- Matrix::rowMeans(px)
    sigma <- sqrt(pmax(Matrix::rowMeans(px^2) - mu^2, 0))
    names(sigma) <- rownames(px)
    scores <- matrix(NA_real_, ncol(x), length(programs))
    coverage <- do.call(rbind, lapply(seq_along(programs), function(j) {
      genes <- unique(as.character(programs[[j]]))
      hit <- intersect(genes, union)
      used <- hit[sigma[hit] > 1e-12]
      ok <- length(used) >= min.program.genes
      if (ok) {
        scores[, j] <<- as.numeric(Matrix::crossprod(
          px[used, , drop = FALSE], 1 / sigma[used]
        )) / length(used) - mean(mu[used] / sigma[used])
      }
      data.frame(
        program = names(programs)[j], requested = length(genes),
        matched = length(hit), variable_genes = length(used),
        status = if (ok) "EVALUABLE" else "NOT_EVALUABLE"
      )
    }))
    if (sum(coverage$status == "EVALUABLE") < 2L) {
      stop("At least two programs need min.program.genes matched, variable genes.")
    }
    assignment <- cellstate_winners(scores, names(programs))
  }
  assignment[is.na(assignment) | !nzchar(assignment)] <- "Unassigned"
  names(assignment) <- colnames(x)
  if (is.null(tree)) {
    if (!is.null(input$embedding)) {
      emb <- input$embedding
      if (is.null(rownames(emb)) || anyDuplicated(rownames(emb)) ||
        !setequal(rownames(emb), colnames(x))) {
        stop("Reduction and expression cell IDs do not match.")
      }
      profiles <- as.matrix(Matrix::crossprod(
        membership, emb[colnames(x), , drop = FALSE]
      )) / sizes
      tree.source <- paste("Pearson distance / complete linkage of cluster means:", reduction)
    } else {
      variance <- rowSums((means - rowMeans(means))^2)
      use <- head(order(variance, decreasing = TRUE), min(500L, sum(variance > 0)))
      profiles <- t(means[use, , drop = FALSE])
      tree.source <- "Pearson distance / complete linkage of up to 500 variable cluster-mean genes"
    }
    rownames(profiles) <- groups
    if (ncol(profiles) < 2L || any(!is.finite(profiles))) {
      stop("Insufficient finite profiles for a tree; supply a labeled hclust.")
    }
    correlation <- stats::cor(t(profiles))
    if (any(!is.finite(correlation))) {
      stop("Cluster correlation is undefined; supply a labeled hclust.")
    }
    tree <- stats::hclust(stats::as.dist(pmax(1 - correlation, 0)), method = "complete")
    tree$labels <- groups
  } else {
    tree.source <- "Caller-supplied hclust"
  }
  if (!inherits(tree, "hclust") || anyDuplicated(tree$labels) ||
    !setequal(tree$labels, groups) || nrow(tree$merge) != length(groups) - 1L ||
    any(!is.finite(tree$height)) || any(tree$height < 0)) {
    stop("tree must be an hclust with unique labels matching all cluster IDs.")
  }
  order <- tree$labels[tree$order]
  clusters <- clusters[match(order, clusters$cluster), , drop = FALSE]
  rownames(clusters) <- NULL
  enrichment <- cellstate_gsea(gsea, groups, pathways, db)
  marker.labels <- cellstate_labels(features, feature.labels, features)
  pathway.labels <- cellstate_labels(
    enrichment$ids, pathway.labels, enrichment$labels
  )
  data <- list(
    clusters = clusters, markers = markers,
    programs = cellstate_composition(meta[[group.by]], assignment),
    samples = cellstate_composition(meta[[group.by]], meta[[sample.by]]),
    pathways = enrichment$data, program.coverage = coverage,
    program.assignment = assignment
  )
  encoded <- cellstate_encode(data, min.pct)
  data$markers <- encoded$markers
  data$pathways <- encoded$pathways
  plot <- cellstate_grob(
    data, tree, features, marker.labels, enrichment$ids, pathway.labels,
    encoded, palette, expression_palette, enrichment_palette,
    palcolor, family, title, panel_label
  )
  provenance <- list(
    input = input$source, input.scale = input.scale, assay = input$assay,
    layer = input$layer, n_cells = ncol(x), n_genes = nrow(x),
    tree = tree.source, gsea = "Precomputed results; no inference run by CellStatePlot"
  )
  if (return_data) {
    return(list(
      plot = plot, data = data, tree = tree, provenance = provenance,
      width_mm = attr(plot, "width_mm"), height_mm = attr(plot, "height_mm")
    ))
  }
  plot
}

cellstate_input <- function(object, assay, layer, reduction) {
  source <- class(object)[1]
  if (is.character(object) && length(object) == 1L) {
    stop("Load the object with readRDS() or h5ad_to_srt() before calling CellStatePlot().")
  }
  if (inherits(object, "Seurat")) {
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
    if (!assay %in% SeuratObject::Assays(object)) stop("Unknown assay: ", assay)
    if (is.null(layer)) layer <- "data"
    if (!layer %in% SeuratObject::Layers(object[[assay]], search = NA)) {
      stop("Select one existing joined layer; join split layers before plotting: ", layer)
    }
    embedding <- if (is.null(reduction)) {
      NULL
    } else {
      if (!reduction %in% SeuratObject::Reductions(object)) {
        stop("Unknown reduction: ", reduction)
      }
      SeuratObject::Embeddings(object[[reduction]])
    }
    return(list(
      x = GetAssayData5(object, assay = assay, layer = layer),
      meta = object[[]], embedding = embedding, gsea.object = object,
      assay = assay, layer = layer, source = source
    ))
  }
  if (!inherits(object, "python.builtin.object") ||
    !reticulate::py_has_attr(object, "obs_names") ||
    !reticulate::py_has_attr(object, "var_names")) {
    stop("object must be a Seurat or Python AnnData object.")
  }
  if (!is.null(assay)) stop("assay applies to Seurat; use layer for AnnData.")
  mat <- if (is.null(layer) || identical(layer, "X")) {
    object$X
  } else {
    if (!layer %in% get_adata_keys(object$layers)) stop("Unknown AnnData layer: ", layer)
    get_adata_element(object$layers, layer)
  }
  if (inherits(mat, "python.builtin.object") &&
    reticulate::py_has_attr(mat, "to_memory")) {
    mat <- mat$to_memory()
  }
  if (inherits(mat, "h5py._hl.dataset.Dataset")) {
    mat <- reticulate::import("numpy", convert = FALSE)$asarray(mat)
  }
  x <- Matrix::t(adata_matrix_to_r(mat))
  dimnames(x) <- list(
    as.character(get_adata_names(object, "var")),
    as.character(get_adata_names(object, "obs"))
  )
  meta <- as.data.frame(py_to_r2(object$obs), check.names = FALSE)
  embedding <- NULL
  if (!is.null(reduction)) {
    if (!reduction %in% get_adata_keys(object$obsm)) stop("Unknown AnnData obsm key: ", reduction)
    embedding <- as.matrix(py_to_r2(get_adata_element(object$obsm, reduction)))
    rownames(embedding) <- colnames(x)
  }
  list(
    x = x, meta = meta, embedding = embedding, gsea.object = NULL,
    assay = NULL, layer = if (is.null(layer)) "X" else layer, source = source
  )
}

cellstate_annotation <- function(meta, group.by, groups, column, map, default) {
  if (!is.null(map)) {
    if (is.null(names(map)) || anyDuplicated(names(map))) {
      stop("Annotation maps need unique cluster-ID names.")
    }
    answer <- as.character(map[groups])
  } else if (!is.null(column)) {
    if (length(column) != 1L) stop("Annotation columns must be single names.")
    answer <- vapply(groups, function(g) {
      value <- unique(as.character(meta[[column]][as.character(meta[[group.by]]) == g]))
      if (length(value) != 1L) stop("Mixed annotations in cluster ", g, "; supply an explicit map.")
      value
    }, character(1))
  } else {
    answer <- rep(default, length.out = length(groups))
  }
  if (anyNA(answer) || any(!nzchar(answer))) stop("Annotations must cover every cluster.")
  unname(answer)
}

cellstate_winners <- function(scores, labels) {
  scores[!is.finite(scores)] <- -Inf
  best <- apply(scores, 1L, max)
  answer <- labels[max.col(scores, ties.method = "first")]
  answer[rowSums(abs(scores - best) <= 1e-12, na.rm = TRUE) > 1L] <- "Ambiguous"
  answer[!is.finite(best)] <- "Unassigned"
  answer
}

cellstate_composition <- function(group, category) {
  tab <- table(as.character(group), as.character(category))
  result <- as.data.frame(tab / rowSums(tab), stringsAsFactors = FALSE)
  names(result) <- c("cluster", "category", "fraction")
  result
}

cellstate_labels <- function(ids, labels, default) {
  if (is.null(labels)) {
    return(unname(default))
  }
  if (is.null(names(labels)) || anyDuplicated(names(labels)) ||
    anyNA(labels[ids])) {
    stop("Label maps must name all displayed IDs uniquely.")
  }
  unname(as.character(labels[ids]))
}

cellstate_gsea <- function(gsea, groups, pathways, db) {
  if (is.list(gsea) && !is.data.frame(gsea)) gsea <- gsea$enrichment
  if (!is.data.frame(gsea)) stop("gsea must be an enrichment table or RunGSEA result.")
  native <- c("Groups", "ID", "NES", "pvalue")
  if (all(native %in% names(gsea))) {
    gsea$cluster <- as.character(gsea$Groups)
    gsea$feature <- as.character(gsea$ID)
    gsea$nes <- gsea$NES
    gsea$p_value <- gsea$pvalue
    if ("p.adjust" %in% names(gsea)) gsea$padj <- gsea$p.adjust
  }
  required <- c("cluster", "feature", "nes", "p_value")
  if (!all(required %in% names(gsea))) {
    stop("gsea needs Groups/ID/NES/pvalue or cluster/feature/nes/p_value.")
  }
  if ("Database" %in% names(gsea)) {
    databases <- unique(as.character(gsea$Database))
    if (is.null(db) && length(databases) > 1L) stop("Select db for a multi-database GSEA table.")
    if (!is.null(db)) {
      if (length(db) != 1L || !db %in% databases) stop("db must select one available GSEA database.")
      gsea <- gsea[as.character(gsea$Database) == db, , drop = FALSE]
    }
  } else if (!is.null(db)) {
    stop("db was supplied but gsea has no Database column.")
  }
  gsea$cluster <- as.character(gsea$cluster)
  gsea$feature <- as.character(gsea$feature)
  if (anyNA(gsea[c("cluster", "feature")]) ||
    any(!nzchar(gsea$cluster)) || any(!nzchar(gsea$feature))) {
    stop("GSEA cluster and pathway IDs must be nonempty.")
  }
  gsea <- gsea[gsea$cluster %in% groups, , drop = FALSE]
  if (!nrow(gsea)) stop("No GSEA cluster IDs match group.by.")
  if (anyDuplicated(gsea[c("cluster", "feature")])) {
    stop("Duplicate GSEA cluster/pathway rows; select one result per combination.")
  }
  if (!is.numeric(gsea$nes) || !is.numeric(gsea$p_value) ||
    any(!is.finite(gsea$nes) & !is.na(gsea$nes)) ||
    any(!is.na(gsea$p_value) & (!is.finite(gsea$p_value) |
      gsea$p_value < 0 | gsea$p_value > 1))) {
    stop("GSEA needs numeric NES and P values between 0 and 1 (or NA).")
  }
  ids <- unique(gsea$feature)
  descriptions <- ids
  if ("Description" %in% names(gsea)) {
    descriptions <- vapply(ids, function(id) {
      values <- unique(as.character(gsea$Description[gsea$feature == id]))
      values <- values[!is.na(values) & nzchar(values)]
      if (length(values) > 1L) stop("Inconsistent GSEA descriptions for ", id)
      if (length(values)) values else id
    }, character(1))
  }
  if (!is.null(pathways)) {
    selected <- vapply(pathways, function(p) {
      if (p %in% ids) {
        return(p)
      }
      hit <- ids[descriptions == p]
      if (length(hit) != 1L) stop("Unknown or ambiguous pathway: ", p)
      hit
    }, character(1))
    ids <- unique(unname(selected))
    descriptions <- descriptions[match(ids, unique(gsea$feature))]
  }
  if (!length(ids)) stop("Select at least one pathway.")
  labels <- gsub("_", " ", sub("^HALLMARK_", "", descriptions))
  labels <- vapply(labels, \(s) paste(strwrap(s, width = 22L), collapse = "\n"), character(1))
  result <- expand.grid(cluster = groups, feature = ids, stringsAsFactors = FALSE)
  index <- match(
    paste(result$cluster, result$feature, sep = "\r"),
    paste(gsea$cluster, gsea$feature, sep = "\r")
  )
  result <- cbind(result, gsea[index, setdiff(names(gsea), c("cluster", "feature")), drop = FALSE])
  rownames(result) <- NULL
  if (!"status" %in% names(result)) result$status <- NA_character_
  result$status <- as.character(result$status)
  valid <- is.finite(result$nes) & is.finite(result$p_value) &
    (is.na(result$status) | result$status != "NOT_EVALUABLE")
  result$status[!valid] <- "NOT_EVALUABLE"
  result$status[valid & is.na(result$status)] <- "EVALUABLE"
  result$neg_log_p <- NA_real_
  result$neg_log_p[valid] <- -log(ifelse(result$p_value[valid] == 0, 1e-10,
    result$p_value[valid]
  ))
  list(data = result, ids = ids, labels = labels)
}
