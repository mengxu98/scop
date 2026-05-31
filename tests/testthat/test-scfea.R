test_that("scFEA bundled species files and annotations are available", {
  human_files <- scop:::scfea_species_files("human")
  mouse_files <- scop:::scfea_species_files("mouse")

  expect_true(all(file.exists(unlist(human_files))))
  expect_true(all(file.exists(unlist(mouse_files))))

  module_info <- scop:::scfea_module_info()
  expect_equal(nrow(module_info), 168)
  expect_true(all(c(
    "module_id", "feature_id", "module_label",
    "reaction_label", "SM_anno"
  ) %in% colnames(module_info)))
  expect_true(all(grepl("^M-", module_info$feature_id)))

  compound_info <- scop:::scfea_compound_info("human")
  expect_equal(nrow(compound_info), 70)
})

test_that("scFEA plotting functions return expected objects", {
  set.seed(1234)
  counts <- matrix(
    rpois(80 * 12, lambda = 5),
    nrow = 80,
    dimnames = list(paste0("Gene", seq_len(80)), paste0("Cell", seq_len(12)))
  )
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- rep(c("Cpa1", "Cre"), each = 6)

  module_info_all <- scop:::scfea_module_info()
  module_info <- module_info_all[
    unique(c(seq_len(8), match("M-150", module_info_all$feature_id))),
    ,
    drop = FALSE
  ]
  flux <- matrix(
    rnorm(nrow(module_info) * 12),
    nrow = nrow(module_info),
    dimnames = list(module_info$feature_id, colnames(srt))
  )
  srt[["scFEAflux"]] <- Seurat::CreateAssayObject(data = flux)
  srt[["scFEAflux"]] <- Seurat::AddMetaData(
    srt[["scFEAflux"]],
    metadata = module_info
  )

  compound_info <- scop:::scfea_compound_info("human")[seq_len(8), , drop = FALSE]
  balance <- matrix(
    rnorm(8 * 12),
    nrow = 8,
    dimnames = list(compound_info$compound_name, colnames(srt))
  )
  srt[["scFEAbalance"]] <- Seurat::CreateAssayObject(data = balance)
  srt[["scFEAbalance"]] <- Seurat::AddMetaData(
    srt[["scFEAbalance"]],
    metadata = compound_info
  )
  srt@tools[["ScFEA"]] <- list(
    flux = flux,
    balance = balance,
    module_info = module_info,
    compound_info = compound_info,
    species = "human",
    n_epoch = 1
  )

  ht <- ScFEAHeatmap(
    srt,
    group.by = "group",
    features = c("M1", "M_2", module_info$reaction_label[3]),
    label_by = "module_reaction",
    add_sm_anno = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    sm_anno_label_rot = 0,
    verbose = FALSE
  )
  expect_s4_class(ht, "Heatmap")
  expect_equal(attr(ht, "width"), 10.8)
  expect_equal(attr(ht, "height"), 11.5)

  expect_equal(
    scop:::scfea_resolve_module_features(
      features = c(
        "M150",
        "M_150",
        "M-150",
        "PRPP -> UMP",
        paste0("M150: PRPP ", intToUtf8(0x2212), "> UMP"),
        paste0("M150: PRPP ", intToUtf8(0x2192), " UMP")
      ),
      module_info = module_info_all,
      available_features = module_info_all$feature_id,
      verbose = FALSE
    ),
    "M-150"
  )

  ht_modules <- ScFEAHeatmap(
    srt,
    group.by = "group",
    modules = "M150: PRPP -> UMP",
    label_by = "module",
    add_sm_anno = FALSE,
    verbose = FALSE
  )
  expect_s4_class(ht_modules, "Heatmap")
  expect_equal(attr(ht_modules, "width"), 7.2)
  expect_equal(attr(ht_modules, "height"), 11.5)

  ht_marked <- ScFEAHeatmap(
    srt,
    group.by = "group",
    modules = module_info$feature_id,
    label_by = "module_reaction",
    mark_features = module_info$feature_id[seq_len(2)],
    heatmap_column_width = grid::unit(10, "mm"),
    column_names_gp = grid::gpar(fontsize = 8),
    column_title_gp = grid::gpar(fontsize = 11),
    row_title_gp = grid::gpar(fontsize = 7.5),
    legend_title_gp = grid::gpar(fontsize = 8),
    legend_labels_gp = grid::gpar(fontsize = 7),
    verbose = FALSE
  )
  expect_s4_class(ht_marked, "Heatmap")

  plots <- ScFEAVolcanoPlot(
    srt,
    group.by = "group",
    ident.1 = "Cpa1",
    ident.2 = "Cre",
    pathways = unique(module_info$SM_anno)[1],
    combine = FALSE,
    verbose = FALSE
  )
  expect_type(plots, "list")
  expect_s3_class(plots[[1]], "ggplot")
  expect_true("data" %in% names(attributes(plots)))
  expect_equal(attr(plots, "width"), 12)
  expect_equal(attr(plots, "height"), 10.4)

  combined_pages <- ScFEAVolcanoPlot(
    srt,
    group.by = "group",
    ident.1 = "Cpa1",
    ident.2 = "Cre",
    pathways = unique(module_info$SM_anno)[1],
    combine = TRUE,
    verbose = FALSE
  )
  expect_type(combined_pages, "list")
  expect_s3_class(combined_pages[[1]], "ggplot")
  expect_silent(ggplot2::ggplot_build(attr(combined_pages, "plots")[[1]]))
  expect_true("data" %in% names(attributes(combined_pages)))
  expect_equal(attr(combined_pages, "width"), 12)
  expect_equal(attr(combined_pages, "height"), 10.4)

  meta_cols_before <- colnames(srt@meta.data)
  auto_pages <- ScFEAVolcanoPlot(
    srt,
    group.by = "group",
    pathways = unique(module_info$SM_anno)[1],
    combine = TRUE,
    verbose = FALSE
  )
  expect_equal(names(auto_pages), c("Cpa1", "Cre"))
  expect_s3_class(auto_pages[["Cpa1"]][[1]], "ggplot")
  expect_silent(ggplot2::ggplot_build(attr(auto_pages[["Cpa1"]], "plots")[[1]]))
  expect_true("data" %in% names(attributes(auto_pages)))
  expect_true(all(c("contrast", "ident.1", "ident.2") %in% colnames(attr(auto_pages, "data"))))
  expect_equal(colnames(srt@meta.data), meta_cols_before)
  expect_error(
    ScFEAVolcanoPlot(
      srt,
      group.by = "group",
      ident.1 = "Cpa1",
      pathways = unique(module_info$SM_anno)[1],
      verbose = FALSE
    ),
    "ident.1"
  )

  balance_plot <- ScFEABalanceBarPlot(
    srt,
    group.by = "group",
    ident.1 = "Cpa1",
    ident.2 = "Cre",
    top_n = 2,
    p_adj_cutoff = 1
  )
  expect_s3_class(balance_plot$plot, "ggplot")
  expect_lte(nrow(balance_plot$data), 4)

  auto_balance <- ScFEABalanceBarPlot(
    srt,
    group.by = "group",
    top_n = 2,
    p_adj_cutoff = 1
  )
  expect_equal(names(auto_balance), c("Cpa1", "Cre"))
  expect_s3_class(auto_balance[["Cpa1"]]$plot, "ggplot")
  expect_true("data" %in% names(attributes(auto_balance)))
  expect_true(all(c("contrast", "ident.1", "ident.2") %in% colnames(attr(auto_balance, "data"))))
  expect_equal(colnames(srt@meta.data), meta_cols_before)
  expect_error(
    ScFEABalanceBarPlot(
      srt,
      group.by = "group",
      ident.1 = "Cpa1",
      top_n = 2,
      p_adj_cutoff = 1
    ),
    "ident.1"
  )
})

test_that("scFEA Python backend imports when dependencies are available", {
  skip_on_os("mac")
  skip_if_not(
    all(vapply(
      c("torch", "numpy", "pandas", "tqdm"),
      reticulate::py_module_available,
      logical(1)
    )),
    "scFEA Python dependencies are not available"
  )

  scfea <- reticulate::import_from_path(
    "scfea",
    path = system.file("python", package = "scop", mustWork = TRUE)
  )
  expect_true("run_scfea" %in% names(scfea))
})

test_that("RunScFEA can run a one-epoch backend smoke test", {
  skip_on_os("mac")
  skip_if_not(
    identical(Sys.getenv("SCOP_RUN_SCFEA_BACKEND_TEST"), "true"),
    "Set SCOP_RUN_SCFEA_BACKEND_TEST=true to run the scFEA training smoke test"
  )
  skip_if_not(
    all(vapply(
      c("torch", "numpy", "pandas", "tqdm"),
      reticulate::py_module_available,
      logical(1)
    )),
    "scFEA Python dependencies are not available"
  )

  genes <- head(scop:::scfea_module_genes("human"), 60)
  counts <- matrix(
    rpois(length(genes) * 6, lambda = 5),
    nrow = length(genes),
    dimnames = list(genes, paste0("Cell", seq_len(6)))
  )
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- RunScFEA(
    srt,
    assay = "RNA",
    layer = "counts",
    species = "human",
    n_epoch = 1,
    verbose = FALSE
  )

  expect_true(all(c("scFEAflux", "scFEAbalance") %in% names(srt@assays)))
  expect_equal(nrow(srt[["scFEAflux"]]), 168)
  expect_equal(nrow(srt[["scFEAbalance"]]), 70)
  expect_identical(colnames(srt[["scFEAflux"]]), colnames(srt))
  expect_true("ScFEA" %in% names(srt@tools))
})
