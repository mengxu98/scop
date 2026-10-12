test_that("SCENICPlus fits TF targets without self predictors", {
  set.seed(82)
  expr <- matrix(rnorm(80 * 4), 4, 80,
                 dimnames = list(c("TF1", "TF2", "G1", "G2"), paste0("c", 1:80)))
  observed <- scenicplus_tf_gene_native(
    expr, c("TF1", "TF2"), targets = "TF1", n_rounds = 40,
    early_stop_window_length = 0, cores = 1, verbose = FALSE
  )
  expected <- grnboost_tree(t(expr), 1:2, 1L, n_rounds = 40,
                            random_seed = 666L, early_stop_window_length = 0L,
                            exclude_self = TRUE)
  expect_equal(observed$importance[observed$TF == "TF2" & observed$target == "TF1"], expected$importance)
  expect_equal(sum(observed$TF == observed$target), 1L)
  expect_false(anyDuplicated(observed[c("TF", "target")]) > 0L)
})

test_that("native GSEA reproduces frozen independent reference statistics", {
  # SCENIC+ 840dab85 / gseapy 0.10.8, 1000 permutations, seed 666.
  observed <- scenicplus_gsea_cpp(list(10:1, 10:1), list(c(1L, 2L, 4L), c(7L, 8L, 10L)), cores = 2)
  expect_equal(observed$NES, c(1.452904068647318, -1.5251544571932922), tolerance = 1e-12)
  expect_equal(observed$pval, c(0.06206896551724138, 0.017361111111111112), tolerance = 1e-12)
  expect_equal(observed$leading_edge, list(c(1L, 2L, 4L), c(7L, 8L, 10L)))
  expect_identical(observed, scenicplus_gsea_cpp(list(10:1, 10:1),
    list(c(1L, 2L, 4L), c(7L, 8L, 10L)), cores = 1))
  empty <- scenicplus_gsea_cpp(list(10:1, 10:1), list(integer(), 1:10))
  expect_true(all(is.na(empty$NES)))
  expect_error(scenicplus_gsea_cpp(list(10:1), list(11L)), "hit")
})

test_that("native ranking ties match frozen PCG reference and ignore scheduling", {
  expr <- cbind(c(3,3,1,0,3,0,2,2,1,0), rep(0,10))
  expected <- rbind(c(1L,2L,5L,8L,0L,7L,4L,3L,6L,9L),
                    c(8L,7L,0L,9L,4L,1L,5L,6L,3L,2L))
  expect_identical(scenicplus_rank_cpp(expr), expected)
  expect_identical(scenicplus_rank_cpp(expr, cores=2), expected)
  expect_error(scenicplus_rank_cpp(replace(expr,1,NA_real_)), "finite")
})

test_that("empty GSEA results cannot silently fall back to unfiltered triplets", {
  tri <- data.frame(TF="TF", region="r",gene="G",score=1)
  tf <- data.frame(TF="TF",target="G",importance=1,rho=1)
  rg <- data.frame(region="r",gene="G",score=1,importance=1,rho=1)
  native <- scenicplus_eregulon_build(tri,tf,rg,min_size=1,min_target_genes=2)
  expect_length(native$regulon_list,0)
  expect_equal(nrow(native$triplets),0L)
  expect_false(native$used_fallback)
  legacy <- scenicplus_eregulon_build(tri,tf,rg,min_size=1,min_target_genes=2,
                                      method="leading_edge")
  expect_length(legacy$regulon_list,1)
  expect_true(legacy$used_fallback)
})

test_that("SCENICPlus self-edge injection is idempotent", {
  expr <- rbind(TF1 = 1:8, TF2 = 8:1, G = c(2, 1, 3, 2, 5, 8, 6, 7))
  adj <- data.frame(TF = c("TF1", "TF1", "TF2"),
                    target = c("G", "TF1", "G"), importance = c(2, 100, 3))
  once <- scenicplus_inject_tf_self(adj, expr)
  twice <- scenicplus_inject_tf_self(once, expr)
  expect_identical(once, twice)
  expect_false(anyDuplicated(once[c("TF", "target")]) > 0L)
  expect_equal(once$importance[once$TF == "TF1" & once$target == "TF1"], 2.00001)
})

test_that("TF model tables retain the frozen reference tie order", {
  regulators <- c("TF3","TF1","TF4","TF2")
  genes <- paste0("G",1:16)
  raw <- data.frame(TF=rep(regulators,16),target=rep(genes,each=4),
    importance=as.numeric((rep(0:15,each=4)+rep(0:3,16))%%3+1),row_id=1:64)
  native <- raw[order(-raw$importance,raw$TF,raw$target),]
  result <- scenicplus_order_tf_edges(native,c(regulators,genes),regulators)
  expect_identical(result$row_id,c(22L,50L,6L,46L,10L,62L,38L,30L,26L,34L,14L,2L,
    54L,58L,18L,42L,32L,52L,56L,16L,4L,44L,64L,28L,8L,12L,60L,40L,24L,36L,20L,48L,
    13L,37L,57L,1L,61L,49L,25L,41L,45L,17L,21L,33L,53L,5L,9L,29L,43L,7L,31L,3L,
    19L,35L,23L,11L,47L,15L,27L,39L,51L,63L,55L,59L))
})

test_that("batch region GBM reuses the scalar model with independent predictors", {
  set.seed(56)
  acc <- matrix(rnorm(75 * 5), 75, 5)
  rna <- cbind(acc[, 1] + acc[, 3], acc[, 2] - acc[, 4], rep(0, 75))
  candidates <- list(c(1L, 3L), c(2L, 4L, 5L), 5L)
  actual <- scenicplus_region_gbm_cpp(acc, rna, candidates, n_rounds = 40, cores = 3)
  serial <- scenicplus_region_gbm_cpp(acc, rna, candidates, n_rounds = 40, cores = 1)
  expect_identical(actual, serial)
  for (g in seq_along(candidates)) {
    regs <- candidates[[g]]
    expr <- cbind(acc[, regs, drop = FALSE], rna[, g])
    ref <- grnboost_tree(expr, seq_along(regs), length(regs) + 1L,
      n_rounds = 40, subsample = 1, random_seed = 666L)
    weights <- numeric(length(regs))
    weights[ref$regulator] <- ref$importance
    if (sum(weights) > 0) weights <- weights / sum(weights)
    expect_equal(actual[[g]], weights, tolerance = 1e-12)
  }
  expect_error(scenicplus_region_gbm_cpp(acc, rna, list(0L, 1L, 2L)), "candidate")
})

test_that("SCENICPlus native fitting never initializes Python", {
  testthat::local_mocked_bindings(.package = "scop",
    run_scenicplus_python = function(...) stop("Python route called"),
    check_python = function(...) stop("Python initialized"),
    scop_python_import = function(...) stop("Python imported"))
  set.seed(7)
  expr <- matrix(rnorm(3 * 40), 3, 40,
                 dimnames = list(c("TF1", "TF2", "G"), paste0("c", 1:40)))
  expect_no_error(scenicplus_tf_gene_native(expr, c("TF1", "TF2"),
    targets = "G", n_rounds = 30, verbose = FALSE))
})

test_that("the native object workflow and its errors stay independent of Python", {
  skip_if_not_installed("Signac")
  blocked <- function(...) stop("Python route must never be called")
  testthat::local_mocked_bindings(.package="scop", run_scenicplus_python=blocked,
    check_python=blocked, scop_python_import=blocked)
  testthat::local_mocked_bindings(.package="reticulate", import=blocked,
    py_run_string=blocked, py_run_file=blocked, py_initialize=blocked,
    initialize_python=blocked)
  original_system2 <- base::system2
  testthat::local_mocked_bindings(.package="base", system2=function(command, ...) {
    if (grepl("python|conda", paste(c(command,unlist(list(...))),collapse=" "),ignore.case=TRUE)) blocked()
    original_system2(command,...)
  })
  cells <- paste0("c",1:30)
  genes <- c("TF",paste0("G",1:39))
  rna <- outer(1:40,1:30,`+`)
  dimnames(rna) <- list(genes,cells)
  regions <- paste0("chr1-",seq(100,2000,100),"-",seq(150,2050,100))
  atac <- outer(1:20,1:30,`*`)
  dimnames(atac) <- list(regions,cells)
  object <- Seurat::CreateSeuratObject(Matrix::Matrix(rna,sparse=TRUE))
  object[["peaks"]] <- Signac::CreateChromatinAssay(Matrix::Matrix(atac,sparse=TRUE),
    sep=c("-","-"),min.cells=0,min.features=0)
  args <- list(object=object,backend="cpp",allow_approximate=TRUE,regulators="TF",
    tf_region_prior=data.frame(TF="TF",region=regions[1:2],score=1),
    region_gene_search_space=expand.grid(region=regions[1:2],gene=genes[2:5]),
    eregulons_prior=data.frame(regulon="TF_direct_+/+_(4g)",target=genes[2:5]),
    min_eregulon_size=2,grn_n_rounds=10,region_gene_n_rounds=10,cores=2,verbose=FALSE)
  fitted <- suppressWarnings(do.call(RunSCENICPlus,args))
  result <- fitted@tools$SCENICPlus
  expect_equal(dim(result$auc),c(30L,1L))
  expect_equal(nrow(result$region_auc),30L)
  expect_false(result$parameters$exact_reference)
  expect_equal(result$parameters$grn_early_stop_window_length,25)
  expect_equal(result$parameters$grn_n_rounds,10)
  args$region_gene_search_space$gene <- "not_in_input"
  expect_error(suppressWarnings(do.call(RunSCENICPlus,args)),"no candidate pairs")
  args$eregulons_prior <- NULL
  args$tf_gene_prior <- data.frame(TF="TF",target=genes[2:7],importance=6:1)
  args$region_gene_prior <- expand.grid(region=regions[1:2],gene=genes[2:5])
  args$region_gene_prior$score <- 1
  args$egrn_min_target_genes <- 2L
  args$egrn_n_perm <- 100L
  from_priors <- suppressWarnings(do.call(RunSCENICPlus,args))@tools$SCENICPlus
  expect_true("rho" %in% names(from_priors$tf_gene))
  expect_true("rho" %in% names(from_priors$region_gene))
  expect_gt(nrow(from_priors$gsea_statistics),0)
})

test_that("group ordering preserves frozen candidate order and duplicates", {
  expect_identical(scenicplus_order_cpp(rep(c("B", "A", "C", "B"), 8)),
    c(2L,30L,6L,26L,10L,22L,14L,18L,1L,29L,28L,25L,24L,21L,20L,17L,
      16L,13L,12L,9L,8L,5L,4L,32L,15L,19L,11L,23L,7L,27L,3L,31L))
  expect_identical(scenicplus_order_cpp(character()), integer())
  expect_error(scenicplus_order_cpp(c("G",NA_character_)), "missing")
})

test_that("region signatures preserve both correlation signs", {
  tri <- data.frame(TF="TF",region=c("r1","r2","r3"),gene=c("G1","G1","G2"),score=1)
  tf <- data.frame(TF="TF", target=c("G1","G2"),rho=c(1,-1))
  rg <- data.frame(region=c("r1","r2","r3"),gene=c("G1","G1","G2"),rho=c(1,-1,1))
  sets <- list("TF_direct_+/+_(1g)"="G1", "TF_direct_+/-_(1g)"="G1",
               "TF_direct_-/+_(1g)"="G2")
  expect_identical(scenicplus_region_signatures(sets,tri,tf,rg),
    list("TF_direct_+/+_(1r)"="r1", "TF_direct_+/-_(1r)"="r2", "TF_direct_-/+_(1r)"="r3"))
})

test_that("tied GSEA rankings preserve the frozen descending permutation", {
  ranking <- stats::setNames(rep(c(4,3,4,2,3,2,1,1),8),paste0("G",1:64))
  expect_setequal(scenicplus_gsea_leading_edge(ranking,paste0("G",1:8)),
    c("G1","G3","G2","G5"))
  expected <- c(1L,17L,59L,57L,51L,49L,43L,41L,35L,27L,25L,19L,33L,9L,3L,11L,
    2L,10L,50L,53L,45L,42L,58L,37L,34L,29L,5L,26L,13L,61L,18L,21L,
    28L,6L,54L,46L,44L,12L,4L,20L,38L,36L,22L,62L,52L,60L,30L,14L,
    56L,55L,63L,39L,48L,47L,7L,8L,40L,15L,32L,31L,24L,23L,16L,64L)
  order <- scenicplus_order_cpp(unname(ranking),decreasing=TRUE)
  expect_identical(order,expected)
  ranking <- ranking[order]
  fit <- scenicplus_gsea_cpp(list(unname(ranking)),list(which(names(ranking)%in%paste0("G",1:8))))
  expect_equal(fit$NES,0.8095510865266078,tolerance=1e-12)
  expect_equal(fit$pval,0.6963226571767497,tolerance=1e-12)
  expect_identical(names(ranking)[fit$leading_edge[[1]]],c("G1","G3","G2","G5","G6","G4"))
  expect_error(scenicplus_order_cpp(c(1,NA_real_)),"finite")
  tf <- data.frame(TF="TF",target=paste0("G",1:64),
    importance=rep(c(4,3,4,2,3,2,1,1),8),rho=1)
  rg <- data.frame(region=paste0("r",1:8),gene=paste0("G",1:8),
    score=1,importance=1,rho=1)
  tri <- data.frame(TF="TF",region=rg$region,gene=rg$gene,score=1)
  assembled <- scenicplus_eregulon_build(tri,tf,rg,min_size=2,min_target_genes=2,
    top_n_region_gene=1L)
  expect_setequal(assembled$regulon_list[[1]],c("G1","G3","G2","G5","G6","G4"))
})

test_that("GBM preserves target precision while casting split features", {
  expr <- cbind(TF=seq_len(50), G=1e9 + seq_len(50) / 10)
  fitted <- grnboost_tree(expr, 1L, 2L, n_rounds=20, subsample=1)
  expect_equal(nrow(fitted), 1L)
  expect_equal(fitted$importance, 1)
})

test_that("constant GBM responses skip trees without changing model edges", {
  expr <- cbind(TF=seq_len(40), variable=seq_len(40)+0.1, zero=0, fixed=1e9)
  fitted <- grnboost_tree_profile(expr,1L,2:4,n_rounds=50,subsample=1)
  reference <- grnboost_tree(expr,1L,2L,n_rounds=50,subsample=1)
  expect_identical(fitted$adjacency,reference)
  expect_equal(unname(fitted$counts[c("targets","rounds","trees")]),c(3L,50L,50L))
  sparse <- grnboost_tree_sparse(methods::as(Matrix::Matrix(expr,sparse=TRUE),"dgCMatrix"),
                                    1L,2:4,n_rounds=50,subsample=1)
  expect_equal(sparse,reference,tolerance=1e-12)
})

test_that("R2G correlations retain ties, constants and missing-value behavior", {
  rna <- rbind(G1=c(1,2,2,4,5), G2=c(1,1,1,1,1))
  acc <- rbind(r1=c(1,1,3,4,5),r2=c(1,NA,2,4,5),r3=5:1)
  hits <- data.frame(region=c("r1","r2","r3","r3"),gene=c("G1","G1","G1","G2"))
  for (mat in list(acc, replace(acc,is.na(acc),2))) {
    result <- scenicplus_rg_gbm(hits,rna,mat,max_region_gene=Inf,n_rounds=15,cores=2)
    expected <- vapply(seq_len(nrow(result)), function(i) {
      suppressWarnings(cor(mat[result$region[i],],rna[result$gene[i],],
        method="spearman",use="complete.obs"))
    },numeric(1))
    expect_equal(result$rho,expected,tolerance=1e-12)
    expect_false("G2" %in% result$gene)
  }
})

test_that("zero-gain parent retains a predictive child split", {
  # Independent Arboreto/sklearn reference: a zero-gain root on TF2,
  # predictive children on TF1, and feature importance c(1, 0).
  expr <- cbind(matrix(rep(c(0,0,0,1,1,0,1,1),10), ncol=2, byrow=TRUE),
                rep(c(0,1,1,0),10))
  for (sparse in c(FALSE, TRUE)) {
    fit <- if (sparse) grnboost_tree_sparse else grnboost_tree
    result <- fit(if (sparse) Matrix::Matrix(expr,sparse=TRUE) else expr, 1:2, 3L,
      n_rounds=1, max_depth=3, max_features=1, subsample=1,
      early_stop_window_length=0, random_seed=666)
    expect_equal(nrow(result), 1L)
    expect_equal(result$regulator, 1L)
    expect_equal(result$importance, 1, tolerance=1e-14)
  }
})
