test_that("tAge preprocessing names are normalized", {
  expect_identical(normalize_tage_model_preprocessing("scaled_diff"), "scaled_diff")
  expect_identical(
    normalize_tage_model_preprocessing(c("scaled_diff", "scaled_diff")),
    "scaled_diff"
  )
  expect_error(
    normalize_tage_model_preprocessing("scaled"),
    "model_preprocessing"
  )
})

test_that("tAge model filenames infer preprocessing names", {
  expect_identical(infer_tage_model_name("EN_Chronoage_scaleddiff.pkl"), "scaled_diff")
  expect_identical(infer_tage_model_name("EN_Chronoage_yugenediff.pkl"), "yugene_diff")
  expect_identical(infer_tage_model_name("EN_Chronoage_yugene.pkl"), "yugene")
})

test_that("fast tAge gene mapping matches upstream mapping", {
  testthat::skip_if_not_installed("tAge")
  testthat::skip_if_not_installed("Biobase")

  metadata_dir <- getFromNamespace("get_metadata_dir", "tAge")()
  gene_table <- utils::read.csv(
    file.path(metadata_dir, "Gene_table_human.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  genes <- head(gene_table[["Gene.Symbol"]][!is.na(gene_table[["Entrez"]])], 120)
  genes <- genes[!duplicated(genes)]
  expr <- matrix(
    seq_len(length(genes) * 4L),
    nrow = length(genes),
    dimnames = list(genes, paste0("sample", seq_len(4L)))
  )
  eset <- Biobase::ExpressionSet(
    assayData = expr,
    phenoData = Biobase::AnnotatedDataFrame(data.frame(
      sample_id = colnames(expr),
      row.names = colnames(expr)
    ))
  )

  expected <- getFromNamespace("map_genes", "tAge")(
    eset,
    species = "human",
    gene_mapping_type = "Gene.Symbol",
    verbose = FALSE
  )
  observed <- tage_map_genes_fast(
    eset,
    species = "human",
    gene_mapping_type = "Gene.Symbol",
    verbose = FALSE
  )

  expect_identical(rownames(observed), rownames(expected))
  expect_equal(Biobase::exprs(observed), Biobase::exprs(expected), tolerance = 0)
  expect_equal(Biobase::fData(observed), Biobase::fData(expected))
})
