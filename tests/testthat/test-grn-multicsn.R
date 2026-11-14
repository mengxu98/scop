grn_multicsn_test_config <- function(mode = "rna", method_args = list()) {
  list(
    method = "multicsn", backend = "r", input_mode = mode,
    assay = "RNA", layer = "data", atac_assay = "ATAC", atac_layer = "data",
    group.by = "celltype", method_args = method_args, resources = list(),
    cores = 1L, seed = 7L, verbose = FALSE
  )
}

grn_multicsn_rna_input <- function() {
  set.seed(19)
  expression <- matrix(stats::rnorm(12L * 48L), 12L,
    dimnames = list(paste0("g", seq_len(12L)), paste0("c", seq_len(48L)))
  )
  expression[5L, ] <- 2 * expression[1L, ] - expression[2L, ] + stats::rnorm(48L, sd = 0.1)
  list(
    object = NULL, expression = expression, group = "all", cells = colnames(expression),
    groups = list(all = colnames(expression)), regulators = rownames(expression)[1:4],
    targets = rownames(expression)[5:8], pseudotime = NULL, input_mode = "rna"
  )
}

grn_multicsn_expect_native <- function(actual, expected) {
  fields <- setdiff(names(expected), "accounting")
  expect_identical(actual[fields], expected[fields])
  expect_identical(names(actual$accounting), names(expected$accounting))
  for (layer in names(expected$accounting)) {
    columns <- setdiff(names(expected$accounting[[layer]]), "seconds")
    expect_identical(actual$accounting[[layer]][columns], expected$accounting[[layer]][columns])
    if ("seconds" %in% names(actual$accounting[[layer]])) {
      expect_true(all(is.finite(actual$accounting[[layer]]$seconds)))
      expect_true(all(actual$accounting[[layer]]$seconds >= 0))
    }
  }
}

grn_multicsn_test_data <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    data <- new.env(parent = emptyenv())
    utils::data("motif2tf", package = "multiCSN", envir = data)
    set.seed(7)
    tfs <- sample(unique(data$motif2tf$tf), 25L)
    targets <- paste0("GENE", seq_len(40L))
    genes <- c(tfs, targets)
    cells <- paste0("cell", seq_len(240L))
    chromosome <- sample(paste0("chr", 1:3), length(genes), TRUE)
    starts <- sample(seq(2e6, 4e7, by = 2e5), length(genes))
    annotation <- GenomicRanges::GRanges(
      chromosome, IRanges::IRanges(starts, width = 5000L),
      strand = sample(c("+", "-"), length(genes), TRUE),
      gene_name = genes, gene_id = paste0("ENSG", seq_along(genes)),
      tx_id = paste0("ENST", seq_along(genes)), gene_biotype = "protein_coding", type = "exon"
    )
    GenomeInfoDb::genome(annotation) <- "hg38"
    peaks <- unique(unlist(lapply(seq_along(genes), function(i) {
      offset <- starts[[i]] + sample(-60000:3000, 4L)
      paste0(chromosome[[i]], "-", offset, "-", offset + 400L)
    })))
    rna <- matrix(stats::rpois(length(genes) * length(cells), 2), length(genes),
      dimnames = list(genes, cells)
    )
    for (j in seq_len(15L)) {
      rna[targets[[j]], ] <- stats::rpois(length(cells), 1 + rna[tfs[[(j %% 5L) + 1L]], ])
    }
    atac <- matrix(stats::rbinom(length(peaks) * length(cells), 3, 0.3), length(peaks),
      dimnames = list(peaks, cells)
    )
    object <- SeuratObject::CreateSeuratObject(
      methods::as(rna, "dgCMatrix"),
      meta.data = data.frame(celltype = rep(c("A", "B"), each = 120L), row.names = cells)
    )
    object[["ATAC"]] <- suppressWarnings(Signac::CreateChromatinAssay(
      methods::as(atac, "dgCMatrix"), ranges = multiCSN::parse_peak_ranges(peaks),
      annotation = annotation, genome = "hg38"
    ))
    object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
    object <- Signac::RunTFIDF(object, assay = "ATAC", verbose = FALSE)
    config <- grn_multicsn_test_config("multiome", list(preset = "state"))
    config$resources$genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
    input <- list(
      object = object,
      expression = SeuratObject::LayerData(object, assay = "RNA", layer = "data"),
      cells = cells, groups = split(cells, object$celltype), group = NULL,
      regulators = NULL, targets = targets, pseudotime = NULL, input_mode = "multiome"
    )
    config$groups <- input$groups
    prepared <- grn_prepare_multicsn(input, config)
    cache <<- list(input = input, prepared = prepared, config = config)
    cache
  }
})

test_that("multiCSN RNA adapter preserves the independent public result", {
  skip_if_not_installed("multiCSN", minimum_version = "0.1.0")
  input <- grn_multicsn_rna_input()
  config <- grn_multicsn_test_config(method_args = list(max_support_size = 3L))
  input <- grn_prepare_multicsn(input, config)
  fit <- infer_multicsn(input, config)
  direct <- multiCSN::inferCSN(
    t(input$expression), regulators = input$regulators, targets = input$targets,
    max_support_size = 3L, cores = 1L, verbose = FALSE
  )
  expect_identical(fit$native, direct)
  expect_identical(fit$cells, input$cells)
  expect_true(any(fit$native$weight < 0))
  expect_identical(fit$layer_specs$tf_gene$sign_column, "weight")
  expect_identical(fit$implementation, "inferCSN::inferCSN")
  input$pseudotime <- matrix(seq(0, 1, length.out = length(input$cells)),
    dimnames = list(input$cells, "trajectory")
  )
  config$method_args$lag_steps <- 1L
  timed <- infer_multicsn(input, config)
  direct_timed <- multiCSN::inferCSN(
    t(input$expression), pseudotime = input$pseudotime,
    regulators = input$regulators, targets = input$targets,
    max_support_size = 3L, lag_steps = 1L, cores = 1L, verbose = FALSE
  )
  expect_identical(timed$native, direct_timed)
  expect_true(timed$parameters$pseudotime)
  config$method_args <- list(preset = "state")
  expect_error(grn_prepare_multicsn(input, config), "Unsupported multiCSN rna")
  input$expression[,] <- 0
  input$pseudotime <- NULL
  config$method_args <- list()
  empty <- infer_multicsn(input, config)
  expect_identical(empty$native, data.frame(regulator = character(), target = character(), weight = numeric()))
  expect_identical(empty$layer_specs$tf_gene$status, "computed")
})

test_that("multiCSN multiome groups preserve complete native endpoints", {
  skip_on_cran()
  skip_if_not_installed("multiCSN", minimum_version = "0.1.0")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  skip_if_not_installed("motifmatchr")
  data <- grn_multicsn_test_data()
  for (group in c("A", "B")) {
    input <- data$prepared
    input$group <- group
    input$cells <- input$groups[[group]]
    input$expression <- input$expression[, input$cells, drop = FALSE]
    fit <- infer_multicsn(input, data$config)
    direct <- multiCSN::fit_layered_network(
      data$prepared$object, cell_group = group, cells = input$cells,
      targets = input$targets, preset = "state", cores = 1L, store = FALSE, verbose = FALSE
    )
    grn_multicsn_expect_native(fit$native, direct)
    expect_identical(fit$cells, input$cells)
    expect_false("object" %in% names(fit$native))
    expect_true(all(c("triplets", "mediated_tf_gene", "accounting", "summary", "settings") %in% names(fit$native)))
    expect_true(all(vapply(fit$layer_specs[1:3], function(spec) identical(spec$sign_column, "standardized_beta"), logical(1))))
    expect_identical(fit$chains$table_path, "triplets")
  }
  input$group <- "A"
  expect_error(infer_multicsn(input, data$config), "not part of")
})

test_that("multiCSN presets retain explicit overrides without stale candidates", {
  skip_on_cran()
  skip_if_not_installed("multiCSN", minimum_version = "0.1.0")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  skip_if_not_installed("motifmatchr")
  data <- grn_multicsn_test_data()
  input <- data$prepared
  input$group <- "B"
  input$cells <- input$groups$B
  input$expression <- input$expression[, input$cells, drop = FALSE]
  config <- data$config
  config$method_args <- list(preset = "timepoint", min_detected = 5L)
  fit <- infer_multicsn(input, config)
  direct <- multiCSN::fit_layered_network(
    data$prepared$object, cell_group = "B", cells = input$cells,
    targets = input$targets, preset = "timepoint", min_detected = 5L,
    cores = 1L, store = FALSE, verbose = FALSE
  )
  grn_multicsn_expect_native(fit$native, direct)
  expect_identical(fit$parameters$min_detected, 5L)
  expect_identical(fit$parameters$tf_region_scope, "motif_peaks")
  expect_true(fit$parameters$renormalize)
  changed <- data$input
  changed$object <- data$prepared$object
  changed$targets <- changed$targets[1:3]
  rebuilt <- grn_prepare_multicsn(changed, data$config)
  expect_identical(rebuilt$multicsn_preparation$targets, changed$targets)
  expect_setequal(multiCSN::get_attribute(rebuilt$object, celltypes = "A", attribute = "genes"), changed$targets)
  expect_true(rebuilt$multicsn_preparation$rebuilt)
  config$layer <- "counts"
  expect_error(grn_prepare_multicsn(data$input, config), "'data' layers")
  config <- data$config
  config$method_args$store <- TRUE
  expect_error(grn_prepare_multicsn(data$input, config), "Unsupported")
  requested <- data$input
  requested$regulators <- rownames(requested$expression)
  filtered <- grn_prepare_multicsn(requested, data$config)
  expect_setequal(filtered$regulators, multiCSN::get_attribute(filtered$object, attribute = "tfs"))
})

test_that("multiCSN layer views preserve independent signs and chain evidence", {
  skip_on_cran()
  skip_if_not_installed("multiCSN", minimum_version = "0.1.0")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  skip_if_not_installed("motifmatchr")
  data <- grn_multicsn_test_data()
  input <- data$prepared
  input$group <- "A"
  input$cells <- input$groups$A
  input$expression <- input$expression[, input$cells, drop = FALSE]
  fit <- infer_multicsn(input, data$config)
  result <- structure(list(schema_version = "1.0", groups = list(A = fit)), class = c("GRNResult", "list"))
  for (layer in c("tf_gene", "tf_region", "region_gene")) {
    edges <- GetGRN(result, group = "A", layer = layer)
    expect_identical(edges$raw_score, fit$native[[layer]]$weight)
    expect_identical(edges$sign, as.integer(sign(fit$native[[layer]]$standardized_beta)))
    expect_identical(edges$standardized_beta, fit$native[[layer]]$standardized_beta)
    expect_identical(edges$deletion_delta_bic, fit$native[[layer]]$deletion_delta_bic)
  }
  expect_identical(GetGRN(result, group = "A", what = "chains"), fit$native$triplets)
  expect_identical(GetGRN(result, group = "A", what = "native"), fit$native)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(result, path)
  expect_identical(GetGRN(readRDS(path), what = "chains"), fit$native$triplets)
})
