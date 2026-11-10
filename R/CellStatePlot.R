#' @title Composite cell-state plot
#'
#' @description
#' Combine a cluster similarity tree, annotations, cell compositions, marker
#' expression and precomputed GSEA in one aligned plot.
#'
#' @md
#' @param object A Seurat or Python AnnData object.
#' @param group.by,sample.by Metadata columns identifying clusters and samples.
#'   Both must be complete. Bars show within-cluster cell fractions.
#' @param features Marker genes in display order. Color shows population-SD
#'   z scores of cluster means, clipped to 0--2; size shows expressing fraction.
#'   Constant genes have z = 0; unmatched genes remain missing.
#' @param program.by Metadata column containing per-cell program assignments.
#' @param program.scores Numeric metadata columns containing program scores.
#'   Optional names label programs. The highest score wins; ties are Ambiguous
#'   and cells without finite scores are Unassigned.
#' @param programs Named gene sets. Each gene is standardized across cells using
#'   population SD; the highest set mean determines assignment. At least two
#'   sets must meet min.program.genes.
#'   Supply exactly one of program.by, program.scores and programs.
#' @param gsea A RunGSEA result or table with Groups, ID, NES and pvalue columns.
#'   Tidy names cluster, feature, nes and p_value are also accepted. NULL uses
#'   stored results for group.by and test.use. Color shows NES clipped to
#'   -2.5--2.5; dot area shows -ln(raw P). Missing or NOT_EVALUABLE entries
#'   use gray crosses.
#' @param pathways Pathway IDs or unambiguous descriptions in display order.
#'   NULL shows all supplied pathways without significance filtering.
#' @param db Value to select from the GSEA table's Database column.
#' @param test.use Test used to locate stored RunGSEA results.
#' @param class.by,subclass.by Metadata columns with one annotation per cluster.
#'   Mixed annotations require an explicit map.
#' @param class.map,subclass.map Named cluster-to-annotation vectors, overriding
#'   the corresponding metadata columns.
#' @param assay Seurat assay. NULL uses DefaultAssay(object).
#' @param layer Expression layer, defaulting to data for Seurat or X for AnnData.
#'   Select one joined Seurat layer.
#' @param input.scale Nonnegative log1p expression or counts. Counts are normalized
#'   to 10,000 per cell and transformed with log1p.
#' @param reduction Seurat reduction or AnnData obsm key used for the tree.
#'   NULL uses up to 500 variable cluster-average gene profiles.
#' @param tree Optional hclust with cluster-ID labels. Otherwise use Pearson
#'   distance and complete linkage on cluster-average profiles.
#' @param min.pct Minimum expressing fraction (0--1) for displaying marker dots.
#' @param expression.cutoff Nonnegative expression threshold used to count
#'   expressing cells, after any counts normalization.
#' @param min.program.genes Minimum matched variable genes per program.
#' @param feature.labels,pathway.labels Named ID-to-label vectors.
#'   Pathway labels support newlines.
#' @param palette Categorical palette passed to thisplot::palette_colors().
#' @param expression_palette,enrichment_palette Continuous palettes for expression
#'   and NES. RdBu maps negative NES to blue and positive NES to red.
#' @param palcolor List of named class, subclass, program and sample colors,
#'   or expression and enrichment color ramps.
#' @param family Font family.
#' @param title,panel_label Optional title and panel letter.
#' @param return_data Include numerical summaries and metadata with the plot.
#' @param flip Arrange clusters in rows with the tree on the left and upright legends.
#'
#' @return A patchwork plot, or a list with the plot and data if return_data = TRUE.
#'
#' @seealso [adata_to_srt()], [h5ad_to_srt()], [RunGSEA()]
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' counts <- GetAssayData5(pancreas_sub, assay = "RNA", layer = "counts")
#' normalized <- counts %*% Matrix::Diagonal(x = 1e4 / Matrix::colSums(counts))
#' dimnames(normalized) <- dimnames(counts)
#' groups <- sort(unique(as.character(pancreas_sub$SubCellType)))
#' scores <- vapply(groups, function(group) {
#'   cells <- colnames(counts)[pancreas_sub$SubCellType == group]
#'   FoldChange(
#'     normalized,
#'     cells.1 = cells, cells.2 = setdiff(colnames(counts), cells),
#'     mean.fxn = function(x) log2(Matrix::rowMeans(x) + 1), fc.name = "log2FC"
#'   )$log2FC
#' }, numeric(nrow(counts)))
#' pathways <- c(
#'   "GO:0031016", "GO:0031017", "GO:0031018", "GO:0030073",
#'   "GO:0008286", "GO:0007049", "GO:0006260", "GO:0006614"
#' )
#' database <- PrepareDB(
#'   species = "Mus_musculus", db = "GO_BP", db_IDtypes = "symbol",
#'   convert_species = FALSE, verbose = FALSE
#' )[["Mus_musculus"]][["GO_BP"]]
#' term2gene <- unique(database$TERM2GENE[
#'   database$TERM2GENE$Term %in% pathways, c("Term", "symbol")
#' ])
#' term2name <- database$TERM2NAME[
#'   database$TERM2NAME$Term %in% pathways, c("Term", "Name")
#' ]
#' set.seed(20261109)
#' gsea <- RunGSEA(
#'   geneID = rep(rownames(counts), length(groups)), geneScore = as.numeric(scores),
#'   geneID_groups = rep(groups, each = nrow(counts)),
#'   species = "Mus_musculus", IDtype = "symbol", result_IDtype = "symbol",
#'   TERM2GENE = term2gene, TERM2NAME = term2name,
#'   minGSSize = 10, maxGSSize = 5000, nPermSimple = 10000, eps = 0,
#'   convert_species = FALSE, verbose = FALSE
#' )
#' CellStatePlot(
#'   pancreas_sub,
#'   group.by = "SubCellType", sample.by = "orig.ident",
#'   features = c("Ins1", "Ins2", "Gcg", "Sst", "Ghrl", "Neurog3", "Fev", "Sox9", "Krt19", "Mki67"),
#'   program.by = "Phase", class.by = "CellType", gsea = gsea, pathways = pathways,
#'   layer = "counts", input.scale = "counts", flip = TRUE
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
  return_data = FALSE,
  flip = FALSE
) {
  input.scale <- match.arg(input.scale)
  if (!is.logical(flip) || length(flip) != 1L || is.na(flip)) {
    stop("flip must be TRUE or FALSE.")
  }
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
    thisutils::log_message(
      "Unmatched markers retained as missing: {paste(absent, collapse = ', ')}",
      message_type = "warning", verbose = TRUE
    )
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
  plot <- cellstate_plot(
    data, tree, features, marker.labels, enrichment$ids, pathway.labels,
    encoded, palette, expression_palette, enrichment_palette,
    palcolor, family, title, panel_label, flip
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
    reticulate::py_get_attr(object, "X")
  } else {
    layers <- reticulate::py_get_attr(object, "layers")
    if (!layer %in% get_adata_keys(layers)) stop("Unknown AnnData layer: ", layer)
    get_adata_element(layers, layer)
  }
  if (inherits(mat, "python.builtin.object") &&
    reticulate::py_has_attr(mat, "to_memory")) {
    mat <- reticulate::py_get_attr(mat, "to_memory")()
  }
  if (inherits(mat, "h5py._hl.dataset.Dataset")) {
    numpy <- reticulate::import("numpy", convert = FALSE)
    mat <- reticulate::py_get_attr(numpy, "asarray")(mat)
  }
  x <- Matrix::t(adata_matrix_to_r(mat))
  dimnames(x) <- list(
    as.character(get_adata_names(object, "var")),
    as.character(get_adata_names(object, "obs"))
  )
  meta <- as.data.frame(py_to_r2(reticulate::py_get_attr(object, "obs")), check.names = FALSE)
  embedding <- NULL
  if (!is.null(reduction)) {
    obsm <- reticulate::py_get_attr(object, "obsm")
    if (!reduction %in% get_adata_keys(obsm)) stop("Unknown AnnData obsm key: ", reduction)
    embedding <- as.matrix(py_to_r2(get_adata_element(obsm, reduction)))
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
  result <- expand.grid(
    cluster = groups, feature = ids, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE
  )
  result <- dplyr::left_join(result, gsea, by = c("cluster", "feature"))
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

cellstate_encode <- function(data, min.pct) {
  markers <- data$markers
  valid <- is.finite(markers$mean_z) & is.finite(markers$proportion)
  markers$show <- valid & markers$proportion >= min.pct
  domain <- if (any(markers$show)) range(markers$proportion[markers$show]) else c(0, 1)
  if (!any(markers$show) || diff(domain) == 0) domain <- c(min.pct, 1)
  if (diff(domain) == 0) domain <- c(0, 1)
  marker.radius <- function(p) {
    u <- pmax(0, pmin(1, (p - domain[1]) / diff(domain)))
    sqrt(0.0354 * (20 + 130 * u^1.8))
  }
  pathways <- data$pathways
  upper <- max(c(0, pathways$neg_log_p), na.rm = TRUE)
  if (upper == 0) upper <- 1
  pathway.radius <- function(p) {
    3.259 * sqrt(pmax(0, pmin(1, p / upper)))
  }
  marker.ticks <- seq(domain[1], domain[2], length.out = 3L)
  format_ticks <- function(ticks) {
    for (digits in 2:17) {
      labels <- format(ticks, trim = TRUE, digits = digits)
      if (!anyDuplicated(labels)) break
    }
    labels
  }
  marker.labels <- format_ticks(marker.ticks)
  size.ticks <- c(.2, .6, 1) * upper
  size.labels <- format(round(size.ticks, 1), trim = TRUE)
  if (anyDuplicated(size.labels)) size.labels <- format_ticks(size.ticks)
  markers$radius <- marker.radius(markers$proportion)
  pathways$radius <- pathway.radius(pathways$neg_log_p)
  list(
    markers = markers, pathways = pathways, marker.valid = valid,
    marker.radius = marker.radius, pathway.radius = pathway.radius,
    marker.ticks = marker.ticks, marker.labels = marker.labels,
    size.ticks = size.ticks, size.labels = size.labels
  )
}

cellstate_palettes <- function(data, palette, expression_palette, enrichment_palette, palcolor) {
  custom <- if (is.null(palcolor)) list() else palcolor
  categories <- c("class", "subclass", "program", "sample")
  if (!is.list(custom) || (length(custom) &&
    (is.null(names(custom)) || anyDuplicated(names(custom)) ||
      any(!names(custom) %in% c(categories, "expression", "enrichment"))))) {
    stop("palcolor must be a named list of class, subclass, program, sample, expression or enrichment colors.")
  }
  for (kind in names(custom)) {
    colors <- custom[[kind]]
    if (!is.character(colors) || !length(colors) || anyNA(colors) ||
      (kind %in% categories && (is.null(names(colors)) ||
        any(!nzchar(names(colors))) || anyDuplicated(names(colors))))) {
      stop("Categorical palcolor entries must be uniquely named color vectors; continuous entries must be color ramps.")
    }
    grDevices::col2rgb(colors)
  }
  categorical <- function(ids, kind) {
    colors <- thisplot::palette_colors(as.character(ids), palette = palette, type = "discrete")
    hit <- intersect(names(custom[[kind]]), names(colors))
    colors[hit] <- custom[[kind]][hit]
    colors
  }
  whiten <- function(color, amount) {
    ramp <- thisplot::palette_colors(
      n = 101, palcolor = c(unname(color), "white"),
      type = "continuous"
    )
    unname(ramp[1L + round(100 * amount)])
  }
  clusters <- data$clusters
  classes <- unique(clusters$class[order(clusters$cluster)])
  base <- categorical(classes, "class")
  subcolors <- stats::setNames(
    rep(NA_character_, length(unique(clusters$subclass))),
    unique(clusters$subclass)
  )
  for (class in classes) {
    ids <- unique(clusters$subclass[clusters$class == class])
    shades <- seq(0, .18, length.out = length(ids))
    for (j in seq_along(ids)) {
      if (is.na(subcolors[ids[j]])) subcolors[ids[j]] <- whiten(base[class], shades[j])
    }
  }
  hit <- intersect(names(custom$subclass), names(subcolors))
  subcolors[hit] <- custom$subclass[hit]
  continuous <- function(name, colors) {
    unname(thisplot::palette_colors(
      n = 257L, palette = name, palcolor = colors,
      type = "continuous"
    ))
  }
  list(
    header = vapply(base, whiten, character(1), amount = .5),
    body = vapply(base, whiten, character(1), amount = .87),
    subclass = subcolors,
    program = categorical(unique(as.character(data$programs$category)), "program"),
    sample = categorical(unique(as.character(data$samples$category)), "sample"),
    expression = continuous(expression_palette, custom$expression),
    enrichment = continuous(enrichment_palette, custom$enrichment)
  )
}

cellstate_plot <- function(data, tree, features, feature.labels, pathways,
                           pathway.labels, encoded, palette, expression_palette,
                           enrichment_palette, palcolor, family, title,
                           panel_label, flip) {
  if (length(family) != 1L || is.na(family)) stop("family must be one font family.")
  if (length(title) > 1L || length(panel_label) > 1L ||
    anyNA(title) || anyNA(panel_label)) {
    stop("title and panel_label must be single strings or NULL.")
  }
  string_width <- if (isTRUE(thisutils::check_r("systemfonts", install = FALSE, verbose = FALSE)[["systemfonts"]])) {
    thisutils::get_namespace_fun("systemfonts", "string_width")
  } else {
    NULL
  }
  measure <- function(s, size = 6, face = "plain") {
    if (!length(s)) {
      return(numeric())
    }
    if (!is.null(string_width)) {
      return(vapply(s, function(value) {
        string_width(
          value,
          family = family, size = size, res = 7200,
          weight = if (face == "bold") "bold" else "normal", italic = face == "italic"
        ) / 100
      }, numeric(1)))
    }
    nchar(s, type = "width") * size * .55
  }
  wrap <- function(s, width) {
    if (grepl("\n", s, fixed = TRUE)) {
      return(s)
    }
    words <- strsplit(s, " ", fixed = TRUE)[[1]]
    lines <- words[1]
    if (length(words) > 1L) {
      for (word in words[-1]) {
        j <- length(lines)
        if (measure(paste(lines[j], word)) <= width) {
          lines[j] <- paste(lines[j], word)
        } else {
          lines <- c(lines, word)
        }
      }
    }
    paste(lines, collapse = "\n")
  }
  palette <- cellstate_palettes(data, palette, expression_palette, enrichment_palette, palcolor)
  n <- nrow(data$clusters)
  xs <- max(28.2732, max(measure(c("Subclasses", "Clusters"))) + 10) +
    (seq_len(n) - 1) * 7.07006
  left <- xs[1] - 5.7712
  right <- tail(xs, 1) + 6.3442
  label.x <- right + 6.0072
  runs <- rle(data$clusters$class)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  edges <- c(left, (xs[head(ends, -1L)] + xs[starts[-1L]]) / 2 - .35, right)
  title.lines <- if (is.null(title) || !nzchar(title)) {
    character()
  } else {
    strsplit(wrap(title, max(110, right - left)), "\n", fixed = TRUE)[[1]]
  }
  display <- function(s) sub("^t cell$", "T cell", gsub("_", " ", s, fixed = TRUE))
  cluster.labels <- display(data$clusters$cluster)
  title.y <- 8
  tree.top <- if (length(title.lines)) title.y + (length(title.lines) - 1L) * 7 + 4 else 4
  tree.height <- 24
  tree.bottom <- tree.top + tree.height
  sub.top <- tree.bottom + 3
  sub.height <- 5
  cluster.top <- sub.top + sub.height + 3
  cluster.bottom <- cluster.top + max(measure(cluster.labels, 5))
  annotation.bottom <- max(cluster.bottom, cluster.top + 5)
  category.layout <- function(colors, top, height, labels = names(colors)) {
    capacity <- max(1L, floor((height - 5) / 6.5) + 1L)
    columns <- min(2L, ceiling(length(colors) / capacity))
    nr <- ceiling(length(colors) / columns)
    column <- (seq_along(colors) - 1L) %/% nr + 1L
    widths <- vapply(seq_len(max(column)), function(j) {
      max(measure(labels[column == j])) + 11.1
    }, numeric(1))
    offsets <- c(0, head(cumsum(widths + 5), -1L))
    list(
      x = label.x + offsets[column],
      y = top + 5 + ((seq_along(colors) - 1L) %% nr) * 6.5,
      right = label.x + sum(widths) + 5 * (length(widths) - 1L),
      bottom = top + 5 + (nr - 1L) * 6.5
    )
  }
  annotation.y <- tree.top + c(4, 11)
  subclass.key <- category.layout(
    palette$subclass, tree.top + 15, annotation.bottom - tree.top - 20,
    display(names(palette$subclass))
  )
  annotation.bottom <- max(annotation.bottom, subclass.key$bottom + 6)
  if (flip) annotation.bottom <- max(cluster.bottom, cluster.top + 5)
  bar.height <- 16
  program.y <- annotation.bottom + 4
  program.key <- category.layout(palette$program, program.y, bar.height)
  sample.y <- if (flip) program.y + bar.height + 4 else max(program.y + bar.height + 4, program.key$bottom + 9)
  sample.key <- category.layout(palette$sample, sample.y, bar.height)
  sample.end <- if (flip) sample.y + bar.height else max(sample.y + bar.height, sample.key$bottom + 2)
  missing <- any(!encoded$marker.valid) || any(!is.finite(data$pathways$neg_log_p))
  gene.top <- sample.end + if (missing) 13 else 4
  gene.span <- (length(features) - 1L) * 5.418856
  gene.height <- max(gene.span + 11.381, max(measure(c("Marker gene", "expression"))) + 4)
  gene.y <- gene.top + (gene.height - gene.span) / 2 +
    (seq_along(features) - 1L) * 5.418856
  gene.bottom <- gene.top + gene.height
  path.lines <- lengths(strsplit(pathway.labels, "\n", fixed = TRUE))
  path.half <- pmax(3.259, path.lines * 2.625)
  path.offsets <- c(0, cumsum(head(path.half, -1L) + tail(path.half, -1L) + 2.5))
  path.content <- tail(path.offsets, 1) + path.half[1] + tail(path.half, 1)
  path.top <- gene.bottom + 3
  path.height <- max(path.content + 4, max(measure(c("Pathway", "enrichment"))) + 4)
  path.y <- path.top + (path.height - path.content) / 2 + path.half[1] + path.offsets
  bottom <- path.top + path.height
  key.offsets <- function(labels) c(0, cumsum(pmax(15, head(measure(labels), -1L) + 8)))
  key.width <- function(label, values) {
    measure(label) + 6 + max(key.offsets(values) + 4 + measure(values))
  }
  key.widths <- c(
    key.width("Expression proportion", encoded$marker.labels),
    key.width("z-normalized expression", c("0", "1", "2")),
    key.width("-ln(P)", encoded$size.labels),
    key.width("GSEA NES", c("-2.5", "0", "2.5"))
  )
  W <- max(
    label.x + 5 + max(measure(unlist(strsplit(c(feature.labels, pathway.labels), "\n", fixed = TRUE)), 5)),
    subclass.key$right + 4, program.key$right + 4, sample.key$right + 4,
    label.x + measure("Cell population") + 4,
    label.x + measure("Program assignment") + 4,
    max(key.widths) + 8
  )
  legend.y <- bottom + 8
  key.columns <- max(key.widths[1:2]) + max(key.widths[3:4]) + 16 <= W
  key.x <- if (key.columns) c(4, 4, 12 + max(key.widths[1:2]), 12 + max(key.widths[1:2])) else rep(4, 4)
  key.y <- legend.y + if (key.columns) c(0, 11, 0, 11) else c(0, 11, 22, 33)
  H <- max(key.y) + 7
  output.width <- W
  output.height <- H
  panel.offset <- 0
  if (flip) {
    output.width <- max(
      bottom + 4, max(key.widths) + 8,
      max(measure(c(display(names(palette$subclass)), names(palette$program), names(palette$sample)))) + 24.1
    )
    if (length(title.lines)) {
      title.lines <- strsplit(wrap(title, output.width - 8), "\n", fixed = TRUE)[[1]]
    }
    panel.offset <- if (length(title.lines) || !is.null(panel_label)) {
      max(10, length(title.lines) * 7 + 4)
    } else {
      0
    }
    landscape.category <- function(colors, top, labels = names(colors),
                                   left = 4, right = output.width - 4) {
      widths <- measure(labels) + 16.1
      xx <- yy <- numeric(length(colors))
      cursor <- left
      row <- top + 9
      for (i in seq_along(colors)) {
        if (cursor > left && cursor + widths[i] > right) {
          cursor <- left
          row <- row + 7
        }
        xx[i] <- cursor
        yy[i] <- row
        cursor <- cursor + widths[i]
      }
      list(x = xx, y = yy, right = max(xx + widths), bottom = max(yy))
    }
    label.height <- max(measure(c(feature.labels, unlist(strsplit(pathway.labels, "\n", fixed = TRUE))), 5))
    legend.y <- label.x + label.height + panel.offset + 10
    annotation.y <- legend.y + c(0, 7)
    subclass.key <- landscape.category(palette$subclass, legend.y + 7, display(names(palette$subclass)))
    program.heading <- subclass.key$bottom + 10
    program.width <- max(sum(measure(names(palette$program)) + 16.1), measure("Program assignment"))
    sample.width <- max(sum(measure(names(palette$sample)) + 16.1), measure("Donors / samples"))
    sample.heading.x <- 4
    if (program.width + sample.width + 12 <= output.width - 8) {
      program.key <- landscape.category(palette$program, program.heading, right = 4 + program.width)
      sample.heading <- program.heading
      sample.heading.x <- 4 + program.width + 12
      sample.key <- landscape.category(palette$sample, sample.heading, left = sample.heading.x)
    } else {
      program.key <- landscape.category(palette$program, program.heading)
      sample.heading <- program.key$bottom + 10
      sample.key <- landscape.category(palette$sample, sample.heading)
    }
    category.bottom <- max(program.key$bottom, sample.key$bottom)
    missing.y <- category.bottom + 10
    key.top <- category.bottom + if (missing) 21 else 12
    key.columns <- max(key.widths[1:2]) + max(key.widths[3:4]) + 16 <= output.width
    key.x <- if (key.columns) c(4, 4, 12 + max(key.widths[1:2]), 12 + max(key.widths[1:2])) else rep(4, 4)
    key.y <- key.top + if (key.columns) c(0, 11, 0, 11) else c(0, 11, 22, 33)
    output.height <- max(key.y) + 7
  }
  colors <- rbind(
    data.frame(
      scale = "expression", value = seq(0, 2, length.out = 257L),
      color = palette$expression
    ),
    data.frame(
      scale = "enrichment", value = seq(-2.5, 2.5, length.out = 257L),
      color = palette$enrichment
    )
  )
  layout <- list(
    W = W, H = H, output.width = output.width, output.height = output.height,
    panel.offset = panel.offset, xs = xs, left = left, right = right, label.x = label.x,
    runs = runs, edges = edges, cluster.labels = cluster.labels,
    title.lines = title.lines, title.y = title.y,
    panel_label = panel_label, family = family,
    tree.top = tree.top, tree.bottom = tree.bottom, tree.height = tree.height,
    sub.top = sub.top, sub.height = sub.height,
    cluster.top = cluster.top, cluster.bottom = cluster.bottom,
    annotation.bottom = annotation.bottom, annotation.y = annotation.y,
    bar.height = bar.height, program.y = program.y, sample.y = sample.y,
    sample.end = sample.end, missing = missing,
    subclass.key = subclass.key, program.key = program.key, sample.key = sample.key,
    gene.y = gene.y, gene.top = gene.top, gene.bottom = gene.bottom,
    path.y = path.y, path.top = path.top, bottom = bottom, legend.y = legend.y,
    key.x = key.x, key.y = key.y, key.offsets = key.offsets,
    measure = measure, display = display
  )
  if (flip) {
    layout$program.heading <- program.heading
    layout$sample.heading <- sample.heading
    layout$sample.heading.x <- sample.heading.x
    layout$missing.y <- missing.y
  }
  state <- list(
    data = data, tree = tree, features = features, feature.labels = feature.labels,
    pathways = pathways, pathway.labels = pathway.labels, encoded = encoded,
    palette = palette, colors = colors, layout = layout, flip = flip
  )
  recorded <- grid::recordGrob(
    grid::grid.draw(render(state)),
    list(state = state, render = cellstate_render),
    name = "cell-state-plot"
  )
  check_r("patchwork", verbose = FALSE)
  plot.theme <- ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  plot <- patchwork::wrap_elements(full = recorded, clip = FALSE) +
    plot.theme
  plot <- patchwork::wrap_plots(plot) + patchwork::plot_annotation(theme = plot.theme)
  attr(plot, "width_mm") <- output.width / 72 * 25.4
  attr(plot, "height_mm") <- output.height / 72 * 25.4
  plot
}

cellstate_render <- function(x) {
  a <- x$layout
  data <- x$data
  palette <- x$palette
  output.width <- a$output.width
  output.height <- a$output.height
  scale <- min(
    grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE) * 72 / output.width,
    grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE) * 72 / output.height
  )
  children <- list()
  add <- function(grob) children[[length(children) + 1L]] <<- grob
  native <- function(value) grid::unit(value, "native")
  rect <- function(xx, yy, w, h, fill, col = NA, lw = .25, panel = TRUE) {
    if (x$flip && panel) {
      original.x <- xx
      xx <- yy
      yy <- original.x + a$panel.offset
      original.w <- w
      w <- h
      h <- original.w
    }
    add(grid::rectGrob(native(xx), native(yy),
      grid::unit(w / output.width, "npc"), grid::unit(h / output.height, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fill = fill, col = col, lwd = lw * scale / .75)
    ))
  }
  line <- function(xx, yy, col = "#D3D3D3", lw = .259, panel = TRUE) {
    if (x$flip && panel) {
      original.x <- xx
      xx <- yy
      yy <- original.x + a$panel.offset
    }
    add(grid::linesGrob(native(xx), native(yy),
      gp = grid::gpar(col = col, lwd = lw * scale / .75, lineend = "square")
    ))
  }
  text <- function(s, xx, yy, size = 6, face = "plain", just = "left",
                   rot = 0, col = "#231F20", panel = TRUE) {
    if (!length(s)) {
      return(invisible(NULL))
    }
    vertical.just <- "bottom"
    if (x$flip && panel) {
      original.x <- xx
      xx <- yy
      yy <- original.x + a$panel.offset
      if (rot == 90) {
        rot <- 0
        just <- switch(just,
          left = "right",
          right = "left",
          just
        )
      } else {
        rot <- -90
        vertical.just <- "top"
      }
    }
    add(grid::textGrob(s, native(xx), native(yy),
      just = c(just, vertical.just), rot = rot,
      gp = grid::gpar(
        fontfamily = a$family, fontsize = size * scale,
        fontface = face, col = col
      )
    ))
  }
  dots <- function(xx, yy, radius, fill, lw = .25, panel = TRUE) {
    if (!length(xx)) {
      return(invisible(NULL))
    }
    yy <- rep(yy, length.out = length(xx))
    fill <- rep(fill, length.out = length(xx))
    if (x$flip && panel) {
      original.x <- xx
      xx <- yy
      yy <- original.x + a$panel.offset
    }
    keep <- is.finite(radius) & radius > 0
    if (!any(keep)) {
      return(invisible(NULL))
    }
    add(grid::circleGrob(native(xx[keep]), native(yy[keep]), native(radius[keep]),
      gp = grid::gpar(fill = fill[keep], col = "black", lwd = lw * scale / .75)
    ))
  }
  cross <- function(xx, yy, panel = TRUE) {
    for (i in seq_along(xx)) {
      line(xx[i] + c(-.7, .7), yy[i] + c(-.7, .7), "#777777", .35, panel = panel)
      line(xx[i] + c(-.7, .7), yy[i] + c(.7, -.7), "#777777", .35, panel = panel)
    }
  }
  continuous <- function(values, kind) {
    if (!length(values)) {
      return(character())
    }
    p <- x$colors[x$colors$scale == kind, ]
    p <- p[order(p$value), ]
    domain <- range(p$value)
    index <- 1L + round((nrow(p) - 1L) *
      pmax(0, pmin(1, (values - domain[1]) / diff(domain))))
    p$color[index]
  }
  ids <- data$clusters$cluster
  n <- length(ids)
  for (i in seq_along(a$runs$values)) {
    class <- a$runs$values[i]
    rect(
      a$edges[i], a$tree.bottom + 1, diff(a$edges)[i],
      a$bottom - a$tree.bottom - 1, palette$body[class]
    )
    rect(a$edges[i], a$tree.top - 1, diff(a$edges)[i], a$tree.height + 2, palette$header[class])
  }
  text(a$panel_label, 1.7457, a$title.y, 8, "bold", panel = FALSE)
  text(a$title.lines, if (x$flip) output.width / 2 else mean(c(a$left, a$right)), a$title.y +
    (seq_along(a$title.lines) - 1L) * 7, just = "centre", panel = FALSE)
  leaf.x <- stats::setNames(a$xs, ids)
  tree.x <- tree.y <- numeric(n - 1L)
  height <- max(x$tree$height)
  for (i in seq_len(n - 1L)) {
    xy <- lapply(x$tree$merge[i, ], function(k) {
      if (k < 0) {
        c(leaf.x[x$tree$labels[-k]], a$tree.bottom)
      } else {
        c(tree.x[k], tree.y[k])
      }
    })
    yy <- a$tree.bottom -
      if (height > 0) a$tree.height * x$tree$height[i] / height else 0
    line(
      c(xy[[1]][1], xy[[1]][1], xy[[2]][1], xy[[2]][1]),
      c(xy[[1]][2], yy, yy, xy[[2]][2]), "#555555", 1.242
    )
    tree.x[i] <- mean(c(xy[[1]][1], xy[[2]][1]))
    tree.y[i] <- yy
  }
  for (j in seq_len(n)) {
    rect(
      a$xs[j] - 2.2436, a$sub.top, 3.536, a$sub.height,
      palette$subclass[data$clusters$subclass[j]]
    )
    text(a$cluster.labels[j], a$xs[j] + 1.063, a$cluster.top +
      a$measure(a$cluster.labels[j], 5), 5, rot = 90, col = "#58595B")
  }
  text("Classes", 16.7994, mean(c(a$tree.top, a$tree.bottom)) + a$measure("Classes") / 2, rot = 90)
  text("Subclasses", a$left - 3, a$sub.top + a$sub.height, just = "right")
  text("Clusters", a$left - 3, a$annotation.bottom, just = "right")
  stack <- function(table, colors, y) {
    for (j in seq_len(n)) {
      cursor <- y + a$bar.height
      for (category in names(colors)) {
        value <- table$fraction[table$cluster == ids[j] & table$category == category]
        if (!length(value) || value == 0) next
        h <- a$bar.height * value
        cursor <- cursor - h
        rect(a$xs[j] - 1.8405, cursor, 3.532, h, colors[category], "black")
      }
    }
  }
  stack(data$programs, palette$program, a$program.y)
  stack(data$samples, palette$sample, a$sample.y)
  for (xx in a$xs) {
    line(c(xx, xx), c(a$gene.top, a$gene.bottom))
    line(c(xx, xx), c(a$path.top, a$bottom - 2.8724))
  }
  for (yy in c(a$gene.y, a$path.y)) line(c(a$left - 3.7733, a$right + 3.2515), rep(yy, 2))
  markers <- data$markers
  mx <- a$xs[match(markers$cluster, ids)]
  my <- a$gene.y[match(markers$feature, x$features)]
  keep <- markers$show
  dots(mx[keep], my[keep], markers$radius[keep], continuous(markers$mean_z[keep], "expression"))
  cross(mx[!x$encoded$marker.valid], my[!x$encoded$marker.valid])
  gsea <- data$pathways
  gx <- a$xs[match(gsea$cluster, ids)]
  gy <- a$path.y[match(gsea$feature, x$pathways)]
  keep <- is.finite(gsea$neg_log_p) & is.finite(gsea$nes)
  dots(gx[keep], gy[keep], gsea$radius[keep], continuous(gsea$nes[keep], "enrichment"))
  cross(gx[!keep], gy[!keep])
  text(x$feature.labels, a$label.x, a$gene.y + 1.5384, 5, "italic", col = "black")
  for (i in seq_along(x$pathway.labels)) {
    words <- strsplit(x$pathway.labels[i], "\n", fixed = TRUE)[[1]]
    yy <- a$path.y[i] + 1.54 - (length(words) - 1L) * 2.625 +
      (seq_along(words) - 1L) * 5.25
    text(words, a$label.x, yy, 5, col = "black")
  }
  text(c("Marker gene", "expression"), c(9.8267, 16.8287),
    mean(c(a$gene.top, a$gene.bottom)) +
      a$measure(c("Marker gene", "expression")) / 2,
    rot = 90
  )
  text(c("Pathway", "enrichment"), c(9.8267, 16.8287),
    mean(c(a$path.top, a$bottom)) +
      a$measure(c("Pathway", "enrichment")) / 2,
    rot = 90
  )
  text(c("Cell population", "annotations"), if (x$flip) 4 else a$label.x, a$annotation.y, panel = FALSE)
  category.key <- function(colors, positions, labels = names(colors)) {
    for (i in seq_along(colors)) {
      xx <- positions$x[i]
      yy <- positions$y[i]
      rect(xx, yy - 3.3437, 8.65, 3.03, colors[i], panel = FALSE)
      text(labels[i], xx + 11.1, yy, panel = FALSE)
    }
  }
  category.key(palette$subclass, a$subclass.key, a$display(names(palette$subclass)))
  text("Program assignment", if (x$flip) 4 else a$label.x,
    if (x$flip) a$program.heading else a$program.y - 2,
    panel = FALSE
  )
  category.key(palette$program, a$program.key)
  text("Donors / samples", if (x$flip) a$sample.heading.x else a$label.x,
    if (x$flip) a$sample.heading else a$sample.y - 2,
    panel = FALSE
  )
  category.key(palette$sample, a$sample.key)
  if (a$missing) {
    cross(if (x$flip) 5 else a$label.x + 1,
      if (x$flip) a$missing.y - 1.6 else a$sample.end + 5.6,
      panel = FALSE
    )
    text("Not evaluable", if (x$flip) 8 else a$label.x + 4,
      if (x$flip) a$missing.y else a$sample.end + 7.2, 5,
      col = "#666666", panel = FALSE
    )
  }
  dot.key <- function(label, labels, radii, fill, xx, yy) {
    text(label, xx, yy + 1.84, col = "black", panel = FALSE)
    dx <- xx + a$measure(label) + 6 + a$key.offsets(labels)
    dots(dx, rep(yy, 3), radii, fill, panel = FALSE)
    text(labels, dx + 4, yy + 1.84, col = "black", panel = FALSE)
  }
  dot.key(
    "Expression proportion", x$encoded$marker.labels,
    x$encoded$marker.radius(x$encoded$marker.ticks), "#D3D3D3", a$key.x[1], a$key.y[1]
  )
  dot.key(
    "z-normalized expression", c("0", "1", "2"), rep(1.505, 3),
    continuous(c(0, 1, 2), "expression"), a$key.x[2], a$key.y[2]
  )
  dot.key(
    "-ln(P)", x$encoded$size.labels,
    x$encoded$pathway.radius(x$encoded$size.ticks), "#333333",
    a$key.x[3], a$key.y[3]
  )
  dot.key(
    "GSEA NES", c("-2.5", "0", "2.5"), rep(1.129, 3),
    continuous(c(-2.5, 0, 2.5), "enrichment"),
    a$key.x[4], a$key.y[4]
  )
  viewport <- grid::viewport(
    width = grid::unit(output.width * scale, "bigpts"), height = grid::unit(output.height * scale, "bigpts"),
    xscale = c(0, output.width), yscale = c(output.height, 0)
  )
  child <- grid::gTree(children = do.call(grid::gList, children), vp = viewport)
  child
}
