test_that("Assay5 variable features are aligned when feature sets differ", {
  skip_if_not_installed("Seurat")
  data(panc8_sub)

  raw_assay <- Seurat::GetAssay(
    Seurat::FindVariableFeatures(panc8_sub, nfeatures = 10, verbose = FALSE)
  )
  append_assay <- Seurat::GetAssay(
    Seurat::FindVariableFeatures(
      subset(panc8_sub, features = rownames(panc8_sub)[1:1000]),
      nfeatures = 5,
      verbose = FALSE
    )
  )

  synced <- scop:::sync_assay_variable_features(
    assay_raw = raw_assay,
    assay_append = append_assay
  )

  expect_equal(
    SeuratObject::VariableFeatures(synced),
    intersect(
      SeuratObject::VariableFeatures(append_assay),
      rownames(raw_assay)
    )
  )
  expect_equal(nrow(synced@meta.data), nrow(raw_assay@meta.data))
})
