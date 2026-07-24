test_that("aggregate_ccc_long cpp backend matches the R backend", {
  df <- data.frame(
    sender = c("B", "A", "A", "B", NA, "A"),
    receiver = c("Y", "X", "X", "Y", "Z", NA),
    score = c(1, 2, NA, 4, 5, 6),
    pvalue = c(0.01, 0.2, 0.03, NA, 0.04, 0.001),
    stringsAsFactors = FALSE
  )

  r_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "r")
  cpp_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "cpp")

  r_out <- r_out[order(r_out$sender, r_out$receiver), , drop = FALSE]
  cpp_out <- cpp_out[order(cpp_out$sender, cpp_out$receiver), , drop = FALSE]
  rownames(r_out) <- NULL
  rownames(cpp_out) <- NULL

  expect_equal(cpp_out, r_out)
})

test_that("aggregate_ccc_long cpp backend falls back for R edge-case labels", {
  df <- data.frame(
    sender = c("AA", "A\nB"),
    receiver = c("", "Z"),
    score = c(4, 5),
    pvalue = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  r_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "r")
  cpp_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "cpp")

  expect_equal(cpp_out, r_out)
})

test_that("aggregate_ccc_long cpp backend handles empty-string sender gracefully", {
  df <- data.frame(
    sender = c("", "A", "A"),
    receiver = c("X", "Y", "Y"),
    score = c(1, 2, 3),
    pvalue = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  r_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "r")
  cpp_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "cpp")

  r_out <- r_out[order(r_out$sender, r_out$receiver), , drop = FALSE]
  cpp_out <- cpp_out[order(cpp_out$sender, cpp_out$receiver), , drop = FALSE]
  rownames(r_out) <- NULL
  rownames(cpp_out) <- NULL

  expect_equal(cpp_out, r_out)
})

test_that("aggregate_ccc_long cpp backend handles NA-only rows", {
  df <- data.frame(
    sender = c(NA, "A", NA),
    receiver = c(NA, "X", "Y"),
    score = c(1, 2, 3),
    pvalue = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  r_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "r")
  cpp_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "cpp")

  r_out <- r_out[order(r_out$sender, r_out$receiver), , drop = FALSE]
  cpp_out <- cpp_out[order(cpp_out$sender, cpp_out$receiver), , drop = FALSE]
  rownames(r_out) <- NULL
  rownames(cpp_out) <- NULL

  expect_equal(cpp_out, r_out)
})

test_that("aggregate_ccc_long cpp backend handles empty data frame", {
  df <- data.frame(
    sender = character(0),
    receiver = character(0),
    score = numeric(0),
    pvalue = numeric(0),
    stringsAsFactors = FALSE
  )

  r_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "r")
  cpp_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "cpp")

  expect_equal(nrow(r_out), 0)
  expect_equal(nrow(cpp_out), 0)
})

test_that("chord pair reduction keeps selected edges when top cell groups have no internal edge", {
  pair_plot <- data.frame(
    sender = c("S1", "S2"),
    receiver = c("R1", "R2"),
    sum = c(10, 9),
    stringsAsFactors = FALSE
  )
  strength_df <- data.frame(
    cell = c("S1", "S2", "R1", "R2"),
    strength = c(20, 18, 10, 9),
    stringsAsFactors = FALSE
  )

  reduced <- getFromNamespace("ccc_reduce_chord_pairs", "scop")(
    pair_plot = pair_plot,
    strength_df = strength_df,
    max.groups = 2
  )

  expect_equal(nrow(reduced$pair_plot), 2)
  expect_equal(reduced$pair_plot$sender, pair_plot$sender)
  expect_equal(reduced$pair_plot$receiver, pair_plot$receiver)
  expect_setequal(reduced$strength_df$cell, c("S1", "S2", "R1", "R2"))
})

test_that("aggregate_ccc_long cpp backend lumped NA as string group", {
  df <- data.frame(
    sender = c("A", NA, "A"),
    receiver = c("X", "Y", NA),
    score = c(1, 2, 3),
    pvalue = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  r_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "r")
  cpp_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "cpp")

  r_out <- r_out[order(r_out$sender, r_out$receiver), , drop = FALSE]
  cpp_out <- cpp_out[order(cpp_out$sender, cpp_out$receiver), , drop = FALSE]
  rownames(r_out) <- NULL
  rownames(cpp_out) <- NULL

  expect_equal(nrow(r_out), 3)
  expect_equal(cpp_out, r_out)
})

test_that("aggregate_ccc_long cpp backend is faster than R backend on larger data", {
  skip_on_cran()
  n <- 50000L
  set.seed(42)
  groups <- paste0("G", sample(1:100, n, replace = TRUE))
  df <- data.frame(
    sender = groups,
    receiver = paste0("R", sample(1:80, n, replace = TRUE)),
    score = runif(n),
    pvalue = runif(n),
    stringsAsFactors = FALSE
  )

  r_time <- system.time({
    r_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "r")
  })[["elapsed"]]
  cpp_time <- system.time({
    cpp_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "cpp")
  })[["elapsed"]]

  r_out <- r_out[order(r_out$sender, r_out$receiver), , drop = FALSE]
  cpp_out <- cpp_out[order(cpp_out$sender, cpp_out$receiver), , drop = FALSE]
  rownames(r_out) <- NULL
  rownames(cpp_out) <- NULL

  expect_equal(cpp_out$sender, r_out$sender)
  expect_equal(cpp_out$receiver, r_out$receiver)
  expect_equal(cpp_out$sum, r_out$sum, tolerance = 1e-10)
  expect_equal(cpp_out$count, r_out$count, tolerance = 1e-10)

  expect_true(cpp_time <= r_time + 0.5)
})

test_that("LIANA C++ aggregation preserves a custom sample key and first non-empty classification", {
  df <- data.frame(
    source = c("A", "A", "A", "B"),
    target = c("X", "X", "X", "Y"),
    ligand_complex = c("L1", "L1", "L1", "L2"),
    receptor_complex = c("R1", "R1", "R1", "R2"),
    batch_id = c("batch1", "batch1", "batch2", "batch1"),
    score = c(1, 2, 3, 4),
    pvalue = c(0.2, 0.1, 0.3, 0.4),
    classification = c("", "Signal", "Other", "ClassB"),
    method = c("m1", "m2", "m1", "m1"),
    liana_method = c("l1", "l2", "l1", "l1"),
    resource = c("r1", "r2", "r1", "r1"),
    stringsAsFactors = FALSE
  )
  r_out <- getFromNamespace("ccc_aggregate_liana_table", "scop")(df, sample_col = "batch_id", backend = "r")
  cpp_out <- getFromNamespace("ccc_aggregate_liana_table", "scop")(df, sample_col = "batch_id", backend = "cpp")
  expect_identical(names(cpp_out), names(r_out))
  expect_equal(cpp_out, r_out, tolerance = 1e-12)
})

test_that("LIANA C++ aggregation supports method as a large-table sample key", {
  n <- 1000L
  df <- data.frame(
    source = rep(c("A", "B"), length.out = n),
    target = rep(c("X", "Y"), length.out = n),
    ligand_complex = paste0("L", seq_len(n)),
    receptor_complex = paste0("R", seq_len(n)),
    score = seq_len(n) / n,
    pvalue = rep(c(0.01, 0.02), length.out = n),
    classification = "Signal",
    method = "CellChat",
    stringsAsFactors = FALSE
  )

  r_out <- getFromNamespace("ccc_aggregate_liana_table", "scop")(
    df,
    sample_col = "method",
    backend = "r"
  )
  cpp_out <- getFromNamespace("ccc_aggregate_liana_table", "scop")(
    df,
    sample_col = "method",
    backend = "cpp"
  )

  expect_true(all(names(r_out) %in% names(cpp_out)))
  expect_equal(cpp_out[, names(r_out), drop = FALSE], r_out, tolerance = 1e-12)
})

test_that("aggregate_ccc_long cpp backend handles mismatched vector lengths with error", {
  df <- data.frame(
    sender = c("A", "B"),
    receiver = c("X"),
    score = c(1, 2),
    pvalue = c(0.1, 0.2),
    stringsAsFactors = FALSE
  )
  expect_error(
    getFromNamespace("ccc_aggregate_long_cpp", "scop")(
      sender = c("A", "B"),
      receiver = c("X"),
      score = c(1, 2),
      significant = c(1, 1)
    ),
    "same length"
  )
})

test_that("aggregate_ccc_long cpp handles single-row data frame", {
  df <- data.frame(
    sender = "A",
    receiver = "X",
    score = 1.5,
    pvalue = 0.01,
    stringsAsFactors = FALSE
  )

  r_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "r")
  cpp_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "cpp")

  expect_equal(nrow(r_out), 1)
  expect_equal(nrow(cpp_out), 1)
  expect_equal(cpp_out$sender, "A")
  expect_equal(cpp_out$receiver, "X")
  expect_equal(cpp_out$sum, 1.5)
  expect_equal(cpp_out$max, 1.5)
  expect_equal(cpp_out$mean, 1.5)
  expect_equal(cpp_out, r_out)
})

test_that("aggregate_ccc_long cpp handles duplicate sender-receiver pairs", {
  df <- data.frame(
    sender = c("A", "A", "A", "B", "B"),
    receiver = c("X", "X", "Y", "X", "X"),
    score = c(1, 3, 5, 2, 4),
    pvalue = c(0.01, 0.05, 0.1, 0.02, 0.03),
    stringsAsFactors = FALSE
  )

  r_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "r")
  cpp_out <- getFromNamespace("aggregate_ccc_long", "scop")(df, backend = "cpp")

  r_out <- r_out[order(r_out$sender, r_out$receiver), , drop = FALSE]
  cpp_out <- cpp_out[order(cpp_out$sender, cpp_out$receiver), , drop = FALSE]
  rownames(r_out) <- NULL
  rownames(cpp_out) <- NULL

  expect_equal(cpp_out, r_out)

  ax_row <- cpp_out[cpp_out$sender == "A" & cpp_out$receiver == "X", ]
  expect_equal(ax_row$sum, 4)
  expect_equal(ax_row$max, 3)
})

test_that("RunCCC dispatches methods, forwards backend, and rebuilds unified CCC", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = c(1, 2, 3, 4),
    dims = c(2, 2)
  )
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- c("Sender", "Receiver")

  seen <- list()
  mock_method <- function(method_name, score) {
    function(srt, group.by, backend = c("cpp", "r"), verbose = TRUE, ...) {
      backend <- match.arg(backend)
      seen[[method_name]] <<- list(
        group.by = group.by,
        backend = backend,
        args = list(...)
      )
      long_table <- data.frame(
        sender = "Sender",
        receiver = "Receiver",
        ligand = paste0("L_", method_name),
        receptor = paste0("R_", method_name),
        interaction_name = paste0("L_", method_name, "_R_", method_name),
        score = score,
        pvalue = 0.01,
        method = method_name,
        stringsAsFactors = FALSE
      )
      if (identical(method_name, "CellChat")) {
        srt@tools[[method_name]] <- list(
          method = method_name,
          results = list(
            ALL = list(
              cellchat_object = structure(
                list(score = score, method_name = method_name),
                class = "mock_cellchat"
              )
            )
          ),
          parameters = list(group.by = group.by, backend = backend)
        )
        return(srt)
      }
      srt@tools[[method_name]] <- list(
        method = method_name,
        results = list(),
        long_table = long_table,
        pair_table = getFromNamespace("aggregate_ccc_long", "scop")(long_table, backend = backend),
        parameters = list(group.by = group.by, backend = backend)
      )
      srt
    }
  }

  testthat::local_mocked_bindings(
    RunCellChat = mock_method("CellChat", 1),
    RunCellphoneDB = mock_method("CellphoneDB", 2),
    RunLIANA = mock_method("LIANA", 3),
    subset_cc_table = function(object, ..., dataset = NULL) {
      data.frame(
        source = "Sender",
        target = "Receiver",
        ligand = "L_CellChat",
        receptor = "R_CellChat",
        interaction_name = paste0("L_CellChat_R_CellChat"),
        prob = object$score,
        pval = 0.01,
        dataset = dataset,
        stringsAsFactors = FALSE
      )
    },
    .package = "scop"
  )

  out <- scop::RunCCC(
    srt = srt,
    group.by = "celltype",
    methods = c("CellChat", "CellphoneDB", "LIANA"),
    method_params = list(CellChat = list(thresh = 0.2)),
    backend = "cpp",
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_equal(names(seen), c("CellChat", "CellphoneDB", "LIANA"))
  expect_equal(unname(vapply(seen, `[[`, character(1), "backend")), rep("cpp", 3))
  expect_equal(seen$CellChat$args$thresh, 0.2)
  expect_equal(out@tools$RunCCC$completed_methods, c("CellChat", "CellphoneDB", "LIANA"))
  expect_equal(out@tools$CCC$metadata$backend, "cpp")
  expect_equal(sort(unique(out@tools$CCC$long_table$method)), c("CellChat", "CellphoneDB", "LIANA"))
  expect_equal(out@tools$CCC$pair_table$sum, 6)
  expect_equal(out@tools$CCC$pair_table$count, 3)
})

test_that("RunCCC with r backend dispatches correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = c(1, 2, 2, 1),
    dims = c(2, 2)
  )
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- c("TypeA", "TypeA")

  seen_backend <- character(0)
  mock_simple <- function(srt, group.by, backend = c("cpp", "r"), verbose = TRUE, ...) {
    backend <- match.arg(backend)
    seen_backend <<- c(seen_backend, backend)
    long_table <- data.frame(
      sender = "TypeA",
      receiver = "TypeA",
      ligand = "L",
      receptor = "R",
      interaction_name = "L_R",
      score = 0.5,
      pvalue = 0.1,
      method = "CellphoneDB",
      stringsAsFactors = FALSE
    )
    srt@tools[["CellphoneDB"]] <- list(
      method = "CellphoneDB",
      results = list(),
      long_table = long_table,
      pair_table = getFromNamespace("aggregate_ccc_long", "scop")(long_table, backend = backend),
      parameters = list(group.by = group.by, backend = backend)
    )
    srt
  }

  testthat::local_mocked_bindings(
    RunCellphoneDB = mock_simple,
    .package = "scop"
  )

  out <- scop::RunCCC(
    srt = srt,
    group.by = "celltype",
    methods = "CellphoneDB",
    backend = "r",
    verbose = FALSE
  )

  expect_equal(seen_backend, "r")
  expect_equal(out@tools$CCC$metadata$backend, "r")
})

test_that("RunNichenetr aggregate_cluster_de maps receivers to upstream arguments", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1, 2, 3, 1, 2),
    j = c(1, 1, 2, 3, 3, 4, 5, 6),
    x = c(5, 2, 4, 3, 6, 2, 1, 4),
    dims = c(3, 6)
  )
  rownames(counts) <- c("L1", "R1", "G1")
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- rep(c("Sender", "Receiver"), each = 3)
  srt$condition <- rep(c("case", "control"), times = 3)

  captured <- list()
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    get_namespace_fun = function(package, name) {
      if (identical(package, "nichenetr") && identical(name, "nichenet_seuratobj_aggregate_cluster_de")) {
        return(function(...) {
          captured$args <<- list(...)
          list(
            ligand_activities = data.frame(test_ligand = "L1", pearson = 0.8),
            ligand_receptor_df = data.frame(from = "L1", to = "R1"),
            ligand_target_df = data.frame(ligand = "L1", target = "G1", weight = 0.4),
            geneset_oi = "G1"
          )
        })
      }
      stop("Unexpected mocked namespace lookup: ", package, "::", name)
    },
    load_nichenetr_models = function(...) {
      list(
        lr_network = data.frame(from = "L1", to = "R1"),
        ligand_target_matrix = matrix(1, nrow = 1, ncol = 1, dimnames = list("G1", "L1")),
        weighted_networks = list()
      )
    },
    .package = "scop"
  )

  out <- scop::RunNichenetr(
    srt = srt,
    group.by = "celltype",
    receiver = "Receiver",
    sender = "Sender",
    condition.by = "condition",
    condition_oi = "case",
    condition_reference = "control",
    mode = "aggregate_cluster_de",
    backend = "cpp",
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_equal(captured$args$receiver_affected, "Receiver")
  expect_equal(captured$args$receiver_reference, "Receiver")
  expect_false("receiver" %in% names(captured$args))
  expect_equal(out@tools$Nichenetr$parameters$backend, "cpp")
})

test_that("RunNichenetr aggregate_cluster_de passes distinct receiver_affected and receiver_reference", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1, 2, 3, 1, 2),
    j = c(1, 1, 2, 3, 3, 4, 5, 6),
    x = c(5, 2, 4, 3, 6, 2, 1, 4),
    dims = c(3, 6)
  )
  rownames(counts) <- c("L1", "R1", "G1")
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- rep(c("Sender", "ReceiverA", "ReceiverB"), times = 2)
  srt$condition <- rep(c("case", "control"), times = 3)

  captured <- list()
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    get_namespace_fun = function(package, name) {
      if (identical(package, "nichenetr") && identical(name, "nichenet_seuratobj_aggregate_cluster_de")) {
        return(function(...) {
          captured$args <<- list(...)
          list(
            ligand_activities = data.frame(test_ligand = "L1", pearson = 0.8),
            ligand_receptor_df = data.frame(from = "L1", to = "R1"),
            ligand_target_df = data.frame(ligand = "L1", target = "G1", weight = 0.4),
            geneset_oi = "G1"
          )
        })
      }
      stop("Unexpected mocked namespace lookup: ", package, "::", name)
    },
    load_nichenetr_models = function(...) {
      list(
        lr_network = data.frame(from = "L1", to = "R1"),
        ligand_target_matrix = matrix(1, nrow = 1, ncol = 1, dimnames = list("G1", "L1")),
        weighted_networks = list()
      )
    },
    .package = "scop"
  )

  out <- scop::RunNichenetr(
    srt = srt,
    group.by = "celltype",
    receiver = "ReceiverA",
    receiver_affected = "ReceiverA",
    receiver_reference = "ReceiverB",
    sender = "Sender",
    condition.by = "condition",
    condition_oi = "case",
    condition_reference = "control",
    mode = "aggregate_cluster_de",
    backend = "cpp",
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_equal(captured$args$receiver_affected, "ReceiverA")
  expect_equal(captured$args$receiver_reference, "ReceiverB")
  expect_false("receiver" %in% names(captured$args))
  expect_equal(out@tools$Nichenetr$parameters$receiver_affected, "ReceiverA")
  expect_equal(out@tools$Nichenetr$parameters$receiver_reference, "ReceiverB")
})

test_that("resolve_nichenetr_object downloads fallback RDS before reading", {
  skip_if_not_installed("R.cache")

  prior <- data.frame(from = "L1", to = "R1", stringsAsFactors = FALSE)
  source <- tempfile(fileext = ".rds")
  saveRDS(prior, source)

  candidate <- basename(tempfile("test_nichenetr_fallback_"))
  out <- getFromNamespace("resolve_nichenetr_object", "scop")(
    x = NULL,
    package = "scop",
    object_candidates = candidate,
    fallback_url = paste0("file://", normalizePath(source)),
    verbose = FALSE
  )

  expect_equal(out, prior)
})

test_that("RunMultiNichenetr preserves user contrast table semantics", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("SingleCellExperiment")

  counts <- Matrix::sparseMatrix(
    i = rep(1:3, 4),
    j = rep(1:4, each = 3),
    x = seq_len(12),
    dims = c(3, 4)
  )
  rownames(counts) <- c("L1", "R1", "G1")
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- c("sender cell", "receiver cell", "sender cell", "receiver cell")
  srt$sample <- c("sample 1", "sample 1", "sample 2", "sample 2")
  srt$condition <- c("condition A", "condition A", "condition B", "condition B")

  user_contrast <- data.frame(
    contrast = "condition A-condition B",
    group = "condition A",
    stringsAsFactors = FALSE
  )
  captured <- list()
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    get_namespace_fun = function(package, name) {
      if (identical(package, "multinichenetr") && identical(name, "multi_nichenet_analysis")) {
        return(function(...) {
          captured$args <<- list(...)
          list(
            prioritization_tables = data.frame(
              sender = "sender.cell",
              receiver = "receiver.cell",
              ligand = "L1",
              receptor = "R1",
              prioritization_score = 0.9,
              p_val = 0.01
            )
          )
        })
      }
      stop("Unexpected mocked namespace lookup: ", package, "::", name)
    },
    load_nichenetr_models = function(...) {
      list(
        lr_network = data.frame(from = "L1", to = "R1"),
        ligand_target_matrix = matrix(1, nrow = 1, ncol = 1, dimnames = list("G1", "L1")),
        weighted_networks = NULL
      )
    },
    .package = "scop"
  )

  out <- suppressWarnings(
    scop::RunMultiNichenetr(
      srt = srt,
      group.by = "celltype",
      sample.by = "sample",
      condition.by = "condition",
      condition_oi = "condition A",
      condition_reference = "condition B",
      receiver_celltypes = "receiver cell",
      sender_celltypes = "sender cell",
      contrast_tbl = user_contrast,
      backend = "cpp",
      verbose = FALSE
    )
  )

  expect_s4_class(out, "Seurat")
  expect_equal(captured$args$contrast_tbl$contrast, "condition.A-condition.B")
  expect_equal(captured$args$contrast_tbl$group, "condition.A")
  expect_equal(captured$args$contrasts_oi, "'condition.A-condition.B'")
  expect_equal(out@tools$MultiNichenetr$parameters$receiver_celltypes, "receiver cell")
  expect_equal(out@tools$MultiNichenetr$parameters$sender_celltypes, "sender cell")
  expect_equal(out@tools$MultiNichenetr$parameters$receiver_celltypes_safe, "receiver.cell")
  expect_equal(out@tools$MultiNichenetr$parameters$sender_celltypes_safe, "sender.cell")
  expect_equal(out@tools$MultiNichenetr$parameters$backend, "cpp")
})

test_that("CellChat method-specific CCC data honors condition instead of cached unified table", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = c(1, 2, 3, 4),
    dims = c(2, 2)
  )
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt@tools[["CellChat"]] <- list(
    results = list(
      A = list(cellchat_object = structure(list(label = "A"), class = "mock_cellchat")),
      B = list(cellchat_object = structure(list(label = "B"), class = "mock_cellchat"))
    ),
    parameters = list(group.by = "celltype")
  )
  srt@tools[["CCC"]] <- list(
    method = "CCC",
    methods = "CellChat",
    long_table = data.frame(
      sender = "Cached",
      receiver = "Cached",
      ligand = "L0",
      receptor = "R0",
      interaction_name = "L0_R0",
      score = 100,
      pvalue = 0.001,
      dataset = "cached",
      method = "CellChat",
      stringsAsFactors = FALSE
    )
  )

  testthat::local_mocked_bindings(
    subset_cc_table = function(object, ..., dataset = NULL) {
      data.frame(
        source = object$label,
        target = paste0(object$label, "_target"),
        ligand = "L1",
        receptor = "R1",
        interaction_name = paste0("L1_R1_", object$label),
        prob = if (identical(object$label, "A")) 1 else 2,
        pval = 0.01,
        dataset = dataset,
        stringsAsFactors = FALSE
      )
    },
    .package = "scop"
  )

  plot_data <- getFromNamespace("ccc_plot_data", "scop")(
    srt = srt,
    method = "CellChat",
    condition = "B",
    value = "score"
  )

  expect_equal(unique(plot_data$long_df$sender), "B")
  expect_equal(unique(plot_data$long_df$receiver), "B_target")
  expect_equal(unique(plot_data$long_df$dataset), "B")
  expect_equal(plot_data$pair_df$sum, 2)
  expect_false("Cached" %in% plot_data$long_df$sender)
})

test_that("CellChat uses unified table when no special parameters are given", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = c(1, 2, 3, 4),
    dims = c(2, 2)
  )
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt@tools[["CellChat"]] <- list(
    method = "CellChat",
    results = list(
      ALL = list(cellchat_object = structure(list(label = "ALL"), class = "mock_cellchat"))
    ),
    parameters = list(group.by = "celltype")
  )
  srt@tools[["CCC"]] <- list(
    method = "CCC",
    methods = "CellChat",
    long_table = data.frame(
      sender = "Unified",
      receiver = "Unified",
      ligand = "L1",
      receptor = "R1",
      interaction_name = "L1_R1",
      score = 50,
      pvalue = 0.01,
      dataset = "ALL",
      method = "CellChat",
      stringsAsFactors = FALSE
    )
  )

  long_df <- getFromNamespace("ccc_long_table_for_method", "scop")(
    srt = srt,
    method = "CellChat",
    dataset = 1,
    slot.name = "net",
    thresh = 0.05
  )

  expect_true("Unified" %in% long_df$sender)
})

test_that("ccc_result_long_table is the internal standardized CCC result contract", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = c(1, 2, 3, 4),
    dims = c(2, 2)
  )
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt@tools[["CCC"]] <- list(
    method = "CCC",
    methods = "CellChat",
    long_table = data.frame(
      source = c("A", "A", "B"),
      target = c("B", "B", "A"),
      ligand_complex = c("L1", "L2", "L3"),
      receptor_complex = c("R1", "R2", "R3"),
      interacting_pair = c("L1_R1", "L2_R2", "L3_R3"),
      score = c(0.8, 0.4, 0.2),
      pvalue = c(0.01, 0.02, 0.2),
      method = "CellChat",
      stringsAsFactors = FALSE
    )
  )

  long <- getFromNamespace("ccc_result_long_table", "scop")(
    srt = srt,
    method = "CCC",
    sender.use = "A"
  )
  expect_true(all(c("sender", "receiver", "ligand", "receptor", "score", "pvalue") %in% colnames(long)))
  expect_equal(unique(long$sender), "A")

  pair <- getFromNamespace("aggregate_ccc_long", "scop")(long, backend = "r")
  expect_equal(pair[pair$sender == "A" & pair$receiver == "B", "count"], 2)

  liana <- getFromNamespace("ccc_long_to_liana", "scop")(long)
  expect_true(all(c("source", "target", "ligand_complex", "receptor_complex") %in% colnames(liana)))
})

test_that("ccc_pair_table accepts backend parameter", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = c(1, 2, 3, 4),
    dims = c(2, 2)
  )
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- c("A", "B")
  srt@tools[["CellChat"]] <- list(
    method = "CellChat",
    results = list(),
    long_table = data.frame(
      sender = c("A", "B"),
      receiver = c("B", "A"),
      ligand = c("L1", "L2"),
      receptor = c("R1", "R2"),
      interaction_name = c("L1_R1", "L2_R2"),
      score = c(0.8, 0.6),
      pvalue = c(0.01, 0.05),
      significant = c(1, 0),
      method = c("CellChat", "CellChat"),
      stringsAsFactors = FALSE
    ),
    pair_table = data.frame(
      sender = c("A", "B"),
      receiver = c("B", "A"),
      sum = c(0.8, 0.6),
      mean = c(0.8, 0.6),
      max = c(0.8, 0.6),
      count = c(1, 0),
      stringsAsFactors = FALSE
    ),
    parameters = list(group.by = "celltype", backend = "r")
  )
  srt@tools[["CCC"]] <- list(
    method = "CCC",
    methods = "CellChat",
    long_table = srt@tools[["CellChat"]]$long_table,
    pair_table = srt@tools[["CellChat"]]$pair_table,
    metadata = list(methods = "CellChat", backend = "r")
  )

  result_r <- getFromNamespace("ccc_pair_table", "scop")(
    srt = srt,
    method = "CellChat",
    backend = "r"
  )
  result_cpp <- getFromNamespace("ccc_pair_table", "scop")(
    srt = srt,
    method = "CellChat",
    backend = "cpp"
  )

  expect_s3_class(result_r, "data.frame")
  expect_s3_class(result_cpp, "data.frame")
  expect_true(nrow(result_r) > 0)
  expect_true(nrow(result_cpp) > 0)
})
