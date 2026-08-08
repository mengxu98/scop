make_ccc_identity_counts <- function() {
  methods::as(Matrix::Diagonal(2), "generalMatrix")
}

test_that("CCC dispatch resolves integrated producers directly", {
  runner <- getFromNamespace("ccc_method_runner", "scop")
  methods <- c(
    "CellChat", "CellphoneDB", "LIANA", "Nichenetr",
    "MultiNichenetr", "SpatialCellChat", "MDIC3"
  )
  for (method in methods) {
    expect_identical(
      runner(method),
      getFromNamespace(paste0("Run", method), "scop")
    )
  }
  expect_identical(runner("CellPhoneDB"), scop::RunCellphoneDB)
  expect_identical(runner("NicheNet"), scop::RunNichenetr)
  expect_error(runner("UnknownCCC"), "Unsupported CCC method")
})

test_that("RunLIANA keeps legacy positional arguments and accepts custom resources", {
  expect_equal(
    names(formals(scop::RunLIANA))[1:9],
    c(
      "srt", "group.by", "method", "resource", "assay", "min_cells",
      "return_all", "backend", "verbose"
    )
  )
  custom <- data.frame(ligand = "L1", receptor = "R1")
  expect_identical(
    getFromNamespace("liana_resolve_resources", "scop")(
      custom, species = "mouse"
    ),
    custom
  )
})

test_that("LIANA resources are discovered from the optional backend", {
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    liana_get_fun = function(fun, package = "liana") {
      expect_equal(package, "liana")
      expect_equal(fun, "show_resources")
      function() c("Consensus", "CellChatDB", "MouseConsensus", "OmniPath")
    },
    .package = "scop"
  )

  resources <- scop::ListCCCDB(db = "LIANA")
  expect_equal(resources$Resource, c("Consensus", "CellChatDB", "MouseConsensus", "OmniPath"))
  expect_equal(resources$Species, c("human", "human", "mouse", "human"))
  expect_equal(resources$Status[resources$Resource == "OmniPath"], "available")
})

test_that("RunLIANA stores official consensus and keeps legacy tables", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = c(2, 1, 3, 4),
    dims = c(2, 2)
  )
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$celltype <- c("A", "B")

  raw <- list(
    natmi = data.frame(
      source = "A", target = "B", ligand.complex = "L1",
      receptor.complex = "R1", score = 0.6, pvalue = 0.2
    ),
    sca = data.frame(
      source = "A", target = "B", ligand.complex = "L1",
      receptor.complex = "R1", score = 0.8, pvalue = 0.1
    )
  )
  consensus <- data.frame(
    source = "A", target = "B", ligand.complex = "L1",
    receptor.complex = "R1", magnitude_rank = 0.05,
    specificity_rank = 0.02
  )
  calls <- character(0)

  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    liana_get_fun = function(fun, package = "liana") {
      if (identical(package, "SingleCellExperiment")) {
        return(function(...) list(...))
      }
      if (identical(package, "utils")) {
        return(function(...) numeric_version("0.1.14"))
      }
      switch(fun,
        show_resources = function() c("Consensus", "MouseConsensus"),
        liana_wrap = function(...) raw,
        rank_aggregate = function(liana_res, ...) {
          calls <<- c(calls, "rank_aggregate")
          consensus
        },
        stop("Unexpected LIANA function: ", fun)
      )
    },
    .package = "scop"
  )

  out <- scop::RunLIANA(
    srt,
    group.by = "celltype",
    method = c("natmi", "sca"),
    species = "human",
    backend = "r",
    verbose = FALSE
  )

  expect_equal(calls, "rank_aggregate")
  expect_true(all(c(
    "results", "long_table", "liana_table", "pair_table",
    "consensus_by_resource", "consensus_table", "primary_table",
    "primary_pair_table"
  ) %in% names(out@tools$LIANA)))
  expect_equal(out@tools$LIANA$consensus_table$magnitude_rank, 0.05)
  expect_equal(out@tools$LIANA$consensus_table$specificity_rank, 0.02)
  expect_equal(out@tools$LIANA$primary_table$score_type, "liana_consensus_priority")
  expect_equal(out@tools$LIANA$primary_table$pvalue_type, "specificity_rank_not_pvalue")
  expect_true(all(is.na(out@tools$LIANA$primary_table$significant)))
  expect_equal(out@tools$LIANA$parameters$consensus, "rank")
  expect_equal(unique(out@tools$LIANA$long_table$resource), "Consensus")
  expect_equal(unique(out@tools$CCC$long_table$score_type), "liana_consensus_priority")
})

test_that("LIANA consensus failure does not mutate the Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- make_ccc_identity_counts()
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$celltype <- c("A", "B")
  before <- names(srt@tools)

  raw <- list(
    natmi = data.frame(source = "A", target = "B", ligand.complex = "L1", receptor.complex = "R1", score = 1),
    sca = data.frame(source = "A", target = "B", ligand.complex = "L1", receptor.complex = "R1", score = 1)
  )
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    liana_get_fun = function(fun, package = "liana") {
      if (identical(package, "SingleCellExperiment")) {
        return(function(...) list(...))
      }
      switch(fun,
        show_resources = function() c("Consensus", "MouseConsensus"),
        liana_wrap = function(...) raw,
        rank_aggregate = function(...) stop("aggregate failed"),
        stop("Unexpected function")
      )
    },
    .package = "scop"
  )

  expect_error(
    scop::RunLIANA(
      srt,
      group.by = "celltype", method = c("natmi", "sca"),
      backend = "r", verbose = FALSE
    ),
    "aggregate failed"
  )
  expect_identical(names(srt@tools), before)
})

test_that("single-method LIANA auto mode does not invent consensus", {
  result <- getFromNamespace("liana_build_consensus", "scop")(
    res = list(natmi = data.frame(x = 1)),
    method = "natmi",
    resource = "Consensus",
    consensus = "auto",
    verbose = FALSE
  )
  expect_equal(result$mode, "none")
  expect_equal(result$status, "not_applicable_single_method")
  expect_equal(nrow(result$table), 0)

  expect_error(
    getFromNamespace("liana_build_consensus", "scop")(
      res = list(natmi = data.frame(x = 1)),
      method = "natmi",
      resource = "Consensus",
      consensus = "rank",
      verbose = FALSE
    ),
    "at least two"
  )
})

test_that("LIANA aggregates multiple resources independently", {
  raw <- list(
    natmi = list(
      Consensus = data.frame(x = 1),
      CellChatDB = data.frame(x = 2)
    ),
    sca = list(
      Consensus = data.frame(x = 3),
      CellChatDB = data.frame(x = 4)
    )
  )
  seen <- character(0)
  testthat::local_mocked_bindings(
    liana_get_fun = function(fun, package = "liana") {
      expect_equal(fun, "rank_aggregate")
      function(liana_res, resource, ...) {
        seen <<- c(seen, resource)
        data.frame(
          source = "A", target = "B", ligand.complex = paste0("L_", resource),
          receptor.complex = "R", magnitude_rank = 0.1,
          specificity_rank = 0.2
        )
      }
    },
    .package = "scop"
  )
  out <- getFromNamespace("liana_build_consensus", "scop")(
    res = raw,
    method = c("natmi", "sca"),
    resource = c("Consensus", "CellChatDB"),
    consensus = "rank",
    verbose = FALSE
  )
  expect_equal(seen, c("Consensus", "CellChatDB"))
  expect_equal(names(out$by_resource), c("Consensus", "CellChatDB"))
  expect_equal(sort(unique(out$table$resource)), c("CellChatDB", "Consensus"))
})

test_that("cross-method combination separates support and visualization rank semantics", {
  df <- data.frame(
    sender = c("A", "A", "A", "A", "A"),
    receiver = c("B", "B", "B", "B", "B"),
    ligand = c("L1", "L1", "L1", "L2", ""),
    receptor = c("R1", "R1", "R1", "R2", ""),
    interaction_name = c("L1_R1", "L1_R1", "L1_R1", "L2_R2", "MDIC3"),
    score = c(10, 9, 0.2, 5, 4),
    pvalue = c(0.01, 0.02, 0.2, 0.03, 1),
    method = c("CellChat", "CellChat", "LIANA", "CellChat", "MDIC3"),
    stringsAsFactors = FALSE
  )

  support <- getFromNamespace("ccc_combine_methods", "scop")(df, "support")
  shared <- support[support$ligand == "L1", , drop = FALSE]
  expect_equal(shared$support_count, 2)
  expect_equal(shared$support_fraction, 1)
  expect_equal(shared$score_type, "method_support_count")
  expect_false(any(support$method == "MDIC3"))

  ranked <- getFromNamespace("ccc_combine_methods", "scop")(df, "rank")
  expect_true(all(ranked$priority_rank >= 0 & ranked$priority_rank <= 1))
  expect_true(all(ranked$support_type == "scop_visualization_consensus"))
  expect_true(all(is.na(ranked$pvalue)))
  expect_equal(ranked$priority_rank[ranked$ligand == "L1"], 0.125)
  expect_true(all(ranked$priority_score >= 0 & ranked$priority_score <= 1))
})

test_that("LIANA export aggregation keeps resources separate", {
  df <- data.frame(
    sender = c("A", "A"), receiver = c("B", "B"),
    ligand = c("L1", "L1"), receptor = c("R1", "R1"),
    score = c(0.7, 0.8), pvalue = c(0.1, 0.2),
    method = "LIANA", resource = c("Consensus", "CellChatDB")
  )
  out <- getFromNamespace("ccc_long_to_liana", "scop")(
    df, aggregate = TRUE
  )
  expect_equal(nrow(out), 2)
  expect_equal(sort(out$resource), c("CellChatDB", "Consensus"))
})

test_that("unified retrieval preserves backend method labels", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  counts <- make_ccc_identity_counts()
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  long <- data.frame(
    sender = c("A", "A"), receiver = c("B", "B"),
    ligand = c("L1", "L1"), receptor = c("R1", "R1"),
    score = c(1, 0.5), pvalue = c(0.01, 0.02),
    method = c("CellChat", "LIANA")
  )
  srt@tools$CCC <- list(long_table = long, primary_table = long)
  expect_equal(
    getFromNamespace("ccc_semantic_long_table", "scop")(
      srt@tools$CCC$long_table, method = NULL
    )$method,
    c("CellChat", "LIANA")
  )
})

test_that("CellChat cached rows honor the requested threshold", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  counts <- make_ccc_identity_counts()
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt@tools$CellChat <- list(primary_table = data.frame(
    sender = c("A", "A"), receiver = c("B", "B"),
    ligand = c("L1", "L2"), receptor = c("R1", "R2"),
    score = c(1, 1), pvalue = c(0.01, 0.08)
  ))
  out <- getFromNamespace("ccc_bundle_long_table", "scop")(
    srt, method = "CellChat", thresh = 0.05
  )
  expect_equal(out$ligand, "L1")
  expect_true(out$significant)
})

test_that("dependent CellphoneDB evidence is disclosed before combination", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 2))
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  long <- data.frame(
    sender = c("A", "A"), receiver = c("B", "B"),
    ligand = c("L1", "L1"), receptor = c("R1", "R1"),
    interaction_name = c("L1_R1", "L1_R1"), score = c(1, 0.8),
    pvalue = c(0.01, 0.02), method = c("CellphoneDB", "LIANA")
  )
  srt@tools$LIANA <- list(parameters = list(method = c("natmi", "cellphonedb")))
  srt@tools$CCC <- list(
    methods = c("CellphoneDB", "LIANA"), long_table = long,
    metadata = list()
  )

  expect_message(
    getFromNamespace("ccc_prepare_combined_object", "scop")(
      srt, method = "CCC", combine_methods = "support"
    ),
    "not independent"
  )
})

test_that("CCC context filters are applied before cross-method combination", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  counts <- make_ccc_identity_counts()
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt@tools$CCC <- list(
    methods = c("CellChat", "LIANA"),
    long_table = data.frame(
      sender = "A", receiver = "B", ligand = "L1", receptor = "R1",
      score = c(1, 0.5), pvalue = c(0.01, 0.02),
      method = c("CellChat", "LIANA"),
      resource = c("Consensus", "Other")
    ),
    metadata = list()
  )
  filtered <- getFromNamespace("ccc_prepare_filtered_object", "scop")(
    srt, method = "CCC", resource = "Consensus"
  )
  combined <- getFromNamespace("ccc_prepare_combined_object", "scop")(
    filtered, method = "CCC", combine_methods = "support"
  )
  expect_equal(combined@tools$CCC$long_table$support_count, 1)
  expect_equal(combined@tools$CCC$long_table$resource, "Consensus")
})

test_that("MultiNicheNet model loading does not build unused weighted networks", {
  testthat::local_mocked_bindings(
    resolve_nichenetr_object = function(x, object_candidates, ...) {
      if (grepl("ligand_target", object_candidates[[1]])) {
        return(matrix(1, nrow = 1, ncol = 1))
      }
      data.frame(from = "L1", to = "R1")
    },
    get_namespace_fun = function(...) {
      stop("construct_weighted_networks must not be resolved")
    },
    .package = "scop"
  )
  out <- getFromNamespace("load_nichenetr_models", "scop")(
    species = "Mus_musculus", include_weighted = FALSE, verbose = FALSE
  )
  expect_null(out$weighted_networks)
})

test_that("RunCCC preflights extended method requirements", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- make_ccc_identity_counts()
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- c("A", "B")

  expect_error(
    scop::RunCCC(
      srt,
      group.by = "celltype", methods = "Nichenetr",
      method_params = list(), verbose = FALSE
    ),
    "receiver"
  )
  expect_error(
    scop::RunCCC(
      srt,
      group.by = "celltype", methods = "Nichenetr",
      method_params = list(Nichenetr = list(receiver = "B")),
      verbose = FALSE
    ),
    "condition.by.*condition_oi.*condition_reference"
  )
  expect_error(
    scop::RunCCC(
      srt,
      group.by = "celltype", methods = "MultiNichenetr",
      method_params = list(MultiNichenetr = list(sample.by = "sample")),
      verbose = FALSE
    ),
    "condition.by"
  )
  expect_error(
    scop::RunCCC(
      srt,
      group.by = "celltype", methods = "SpatialCellChat",
      method_params = list(), verbose = FALSE
    ),
    "spatial image"
  )
  expect_error(
    scop::RunCCC(
      srt,
      group.by = "celltype", methods = "MDIC3",
      method_params = list(), verbose = FALSE
    ),
    "grn.*grn_method"
  )
})

test_that("RunCCC skip_failed records preflight failures and continues", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  counts <- make_ccc_identity_counts()
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- c("A", "B")

  testthat::local_mocked_bindings(
    RunCellChat = function(srt, ...) {
      srt@tools$CellChat <- list(primary_table = data.frame(
        sender = "A", receiver = "B", ligand = "L1", receptor = "R1",
        score = 1, pvalue = 0.01, method = "CellChat"
      ))
      srt
    },
    .package = "scop"
  )
  out <- suppressWarnings(scop::RunCCC(
    srt,
    group.by = "celltype", methods = c("CellChat", "Nichenetr"),
    skip_failed = TRUE, backend = "r", verbose = FALSE
  ))
  expect_equal(out@tools$RunCCC$status$status, c("completed", "failed"))
  expect_match(out@tools$RunCCC$status$message[2], "receiver")
})

test_that("RunCCC dispatches all four design-specific supported methods", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(1, 1), dims = c(2, 2)
  )
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- c("A", "B")
  srt$sample <- c("S1", "S2")
  srt$condition <- c("case", "control")
  srt$col <- c(1, 2)
  srt$row <- c(1, 2)

  seen <- list()
  mock_srt_method <- function(method) {
    function(srt, group.by, backend = "r", verbose = TRUE, ...) {
      seen[[method]] <<- list(group.by = group.by, backend = backend, args = list(...))
      long <- data.frame(
        sender = "A", receiver = "B", ligand = paste0("L_", method),
        receptor = paste0("R_", method), interaction_name = paste0("L_R_", method),
        score = 1, pvalue = 0.01, method = method, stringsAsFactors = FALSE
      )
      srt@tools[[method]] <- list(
        method = method, long_table = long, primary_table = long,
        pair_table = getFromNamespace("aggregate_ccc_long", "scop")(long, backend = "r"),
        parameters = list(group.by = group.by)
      )
      srt
    }
  }
  mock_mdic3 <- function(object, group.by, verbose = TRUE, ...) {
    seen$MDIC3 <<- list(group.by = group.by, args = list(...))
    long <- data.frame(
      sender = "A", receiver = "B", ligand = "L_MDIC3", receptor = "R_MDIC3",
      interaction_name = "L_R_MDIC3", score = 1, pvalue = NA_real_,
      method = "MDIC3", stringsAsFactors = FALSE
    )
    object@tools$MDIC3 <- list(
      method = "MDIC3", long_table = long, primary_table = long,
      pair_table = getFromNamespace("aggregate_ccc_long", "scop")(long, backend = "r"),
      parameters = list(group.by = group.by)
    )
    object
  }

  testthat::local_mocked_bindings(
    RunNichenetr = mock_srt_method("Nichenetr"),
    RunMultiNichenetr = mock_srt_method("MultiNichenetr"),
    RunSpatialCellChat = mock_srt_method("SpatialCellChat"),
    RunMDIC3 = mock_mdic3,
    .package = "scop"
  )

  out <- scop::RunCCC(
    srt,
    group.by = "celltype",
    methods = c("Nichenetr", "MultiNichenetr", "SpatialCellChat", "MDIC3"),
    method_params = list(
      Nichenetr = list(
        receiver = "B", condition.by = "condition",
        condition_oi = "case", condition_reference = "control"
      ),
      MultiNichenetr = list(
        sample.by = "sample", condition.by = "condition",
        condition_oi = "case", condition_reference = "control",
        receiver_celltypes = "B"
      ),
      SpatialCellChat = list(coord.cols = c("col", "row")),
      MDIC3 = list(grn_method = "correlation")
    ),
    backend = "r",
    verbose = FALSE
  )

  expect_equal(
    names(seen),
    c("Nichenetr", "MultiNichenetr", "SpatialCellChat", "MDIC3")
  )
  expect_false("backend" %in% names(seen$MDIC3$args))
  expect_equal(
    sort(unique(out@tools$CCC$long_table$method)),
    sort(c("Nichenetr", "MultiNichenetr", "SpatialCellChat", "MDIC3"))
  )
})

test_that("unified CCC support aggregation combines stored methods", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 2))
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  long <- data.frame(
    sender = c("A", "A"), receiver = c("B", "B"),
    ligand = c("L1", "L1"), receptor = c("R1", "R1"),
    interaction_name = c("L1_R1", "L1_R1"),
    score = c(2, 0.8), pvalue = c(0.01, 0.02),
    method = c("CellphoneDB", "LIANA"), stringsAsFactors = FALSE
  )
  for (method in c("CellphoneDB", "LIANA")) {
    method_long <- long[long$method == method, , drop = FALSE]
    srt@tools[[method]] <- list(
      method = method, long_table = method_long, primary_table = method_long,
      parameters = list(group.by = "celltype")
    )
  }
  srt@tools$CCC <- list(
    method = "CCC", methods = c("CellphoneDB", "LIANA"),
    long_table = long, metadata = list()
  )

  support <- getFromNamespace("ccc_prepare_combined_object", "scop")(
    srt, method = "CCC", combine_methods = "support"
  )
  expect_equal(unique(support@tools$CCC$long_table$score_type), "method_support_count")
})

test_that("CCC discovery and access adapt v1 long tables without rewriting objects", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 2))
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  legacy <- data.frame(
    sender = "A", receiver = "B", ligand = "L1", receptor = "R1",
    interaction_name = "L1_R1", score = 0.7, pvalue = 0.02,
    stringsAsFactors = FALSE
  )
  srt@tools$CellphoneDB <- list(
    method = "CellphoneDB", long_table = legacy,
    metadata = list()
  )

  before <- srt@tools$CellphoneDB$long_table
  adapted <- getFromNamespace("ccc_semantic_long_table", "scop")(
    srt@tools$CellphoneDB$long_table, method = "CellphoneDB"
  )
  expect_true(all(c(
    "method", "resource", "score_type", "score_direction",
    "priority_rank", "priority_score", "pvalue_type", "support_type",
    "producer", "backend_version"
  ) %in% colnames(adapted)))
  expect_equal(adapted$method, "CellphoneDB")
  expect_identical(srt@tools$CellphoneDB$long_table, before)
})

test_that("stored CellChat primary tables enter unified results without native re-extraction", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 2))
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  primary <- data.frame(
    sender = "A", receiver = "B", ligand = "L1", receptor = "R1",
    interaction_name = "L1_R1", score = 0.8, pvalue = 0.01,
    stringsAsFactors = FALSE
  )
  srt@tools$CellChat <- list(
    results = list(), primary_table = primary, long_table = primary,
    metadata = list()
  )

  unified <- getFromNamespace("ccc_build_unified_bundle", "scop")(
    srt, methods = "CellChat"
  )
  expect_equal(nrow(unified$long_table), 1L)
  expect_equal(unified$long_table$method, "CellChat")
  expect_equal(
    srt@tools$CellChat$primary_table$score,
    0.8
  )
})

test_that("LIANA result accessor preserves official consensus ranks", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 2))
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  consensus <- data.frame(
    sender = "A", receiver = "B", ligand = "L1", receptor = "R1",
    interaction_name = "L1_R1", magnitude_rank = 0.05,
    specificity_rank = 0.02, resource = "Consensus"
  )
  srt@tools$LIANA <- list(
    method = "LIANA", consensus_table = consensus,
    consensus_by_resource = list(Consensus = consensus),
    primary_table = consensus, long_table = consensus
  )

  expect_equal(
    srt@tools$LIANA$consensus_by_resource$Consensus$magnitude_rank,
    0.05
  )
})

test_that("native access never substitutes a raw result", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 2))
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt@tools$CellphoneDB <- list(
    long_table = data.frame(
      sender = "A", receiver = "B", ligand = "L1", receptor = "R1",
      score = 1, pvalue = 0.01
    ),
    raw_result = list(not_native = TRUE)
  )

  expect_null(
    srt@tools$CellphoneDB$native_object %||%
      srt@tools$CellphoneDB$cellchat_object
  )
  expect_true(srt@tools$CellphoneDB$raw_result$not_native)

  srt@tools$MDIC3 <- list(
    mdic3_matrix = matrix(1, 1, 1),
    cellular_communication = matrix(2, 1, 1),
    cellular_communication_log = matrix(3, 1, 1),
    celltype_communication_raw = matrix(4, 1, 1),
    celltype_communication = matrix(5, 1, 1),
    long_table = data.frame(
      sender = "A", receiver = "B", ligand = "", receptor = "",
      score = 5, pvalue = 1
    )
  )
  expect_named(
    srt@tools$MDIC3[intersect(
      c(
        "mdic3_matrix", "cellular_communication",
        "cellular_communication_log", "celltype_communication_raw",
        "celltype_communication"
      ),
      names(srt@tools$MDIC3)
    )],
    c(
      "mdic3_matrix", "cellular_communication",
      "cellular_communication_log", "celltype_communication_raw",
      "celltype_communication"
    )
  )

  srt@tools$SpatialCellChat <- list(long_table = data.frame(
    sender = "A", receiver = "B", ligand = "L1", receptor = "R1",
    score = 0.5, pvalue = 0.01
  ))
  expect_equal(
    nrow(getFromNamespace("aggregate_ccc_long", "scop")(
      srt@tools$SpatialCellChat$long_table, backend = "r"
    )),
    1L
  )
})

test_that("unified resource and sample filters are exact and explicit", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 2))
  rownames(counts) <- c("L1", "R1")
  colnames(counts) <- c("Cell1", "Cell2")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  long <- data.frame(
    sender = "A", receiver = "B", ligand = c("L1", "L2"),
    receptor = c("R1", "R2"), interaction_name = c("L1_R1", "L2_R2"),
    score = c(0.8, 0.6), pvalue = c(0.01, 0.02), method = "LIANA",
    resource = c("Consensus", "CellChatDB"), sample = c("S1", "S2")
  )
  srt@tools$LIANA <- list(long_table = long, primary_table = long)
  srt@tools$CCC <- list(
    methods = "LIANA", long_table = long,
    metadata = list()
  )

  filtered <- getFromNamespace("ccc_prepare_filtered_object", "scop")(
    srt, method = "LIANA", resource = "Consensus", sample = "S1"
  )
  expect_equal(nrow(filtered@tools$CCC$long_table), 1)
  expect_equal(filtered@tools$CCC$long_table$ligand, "L1")

  no_provenance <- srt
  no_provenance@tools$CCC$long_table$resource <- NA_character_
  expect_error(
    getFromNamespace("ccc_prepare_filtered_object", "scop")(
      no_provenance, method = "CCC", resource = "Consensus"
    ),
    "does not provide.*resource"
  )
})
