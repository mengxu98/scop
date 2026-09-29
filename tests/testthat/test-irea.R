make_irea_fixture <- function() {
  mat <- matrix(0, nrow=4, ncol=12,
                dimnames=list(c("G1","G2","G3","G4"),paste0("cell",seq_len(12))))
  mat["G1",] <- c(rep(1,4),rep(30,4),rep(2,4))
  mat["G2",] <- c(rep(10,4),rep(1,4),rep(10,4))
  mat["G3",] <- 2
  mat["G4",] <- c(rep(2,8),rep(30,4))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(mat,sparse=TRUE))
  object <- Seurat::NormalizeData(object,verbose=FALSE)
  object$sample <- rep(c("PBS","IL12","IL15"),each=4)
  object$polarization <- c(rep("None",8),rep("NK-e",4))
  structure(list(object=object,
    cytokine=data.frame(Cytokine_Str=c("IL12","IL15"),Gene=c("G1","G4"),FDR=c(0.001,0.001)),
    polarization=data.frame(Polarization="NK-e",Gene="G4",P_adj=0.001),
    cell_type="NK_cell",species="Mouse",paths=character(),checksum_md5=character(),
    provenance="synthetic test"),class="irea_reference")
}

test_that("gene scores and projection keep their response direction", {
  r <- make_irea_fixture()
  a <- RunIREA(r,genes="G1")$table
  expect_gt(a$effect[a$term=="IL12"],0)
  expect_lt(a$effect[a$term=="IL15"],a$effect[a$term=="IL12"])
  p <- RunIREA(r,matrix=c(G1=1,G2=-1),gene_diff_cutoff=0)$table
  expect_gt(p$effect[p$term=="IL12"],0)
  n <- RunIREA(r,matrix=c(G1=-1,G2=1),gene_diff_cutoff=0)$table
  expect_lt(n$effect[n$term=="IL12"],0)
  expect_equal(p$fdr,stats::p.adjust(p$p_value,"BH"))
})

test_that("hypergeometric and polarization retain separate meanings", {
  r <- make_irea_fixture()
  h <- RunIREA(r,genes="G1",method="hypergeometric")$table
  expect_equal(h$effect[h$term=="IL12"],1)
  expect_equal(h$effect[h$term=="IL15"],0)
  s <- RunIREA(r,genes="G4",analysis="cell_polarization")$table
  expect_gt(s$effect,0)
  expect_s3_class(IREAPlot(s <- RunIREA(r,genes="G4",analysis="cell_polarization"),"radar"),"ggplot")
})

test_that("invalid input fails instead of returning a biological zero", {
  r <- make_irea_fixture()
  expect_error(RunIREA(r,genes="UNKNOWN"),"No input genes match")
  expect_error(RunIREA(r,matrix=c(G1=0,G2=0)),"nonzero")
  expect_error(RunIREA(r,genes="G1",matrix=c(G1=1)),"exactly one")
  expect_error(RunIREA(r,matrix=c(G1=1,G1=2)),"duplicate")
})

test_that("Seurat contrasts use the same projection core", {
  r <- make_irea_fixture()
  o <- r$object
  out <- RunIREA(r,object=o,group_by="sample",case="IL12",control="PBS",
                 gene_diff_cutoff=0)
  expect_s4_class(out,"Seurat")
  expect_s3_class(out@tools$IREA,"irea_result")
  expect_true(all(is.finite(out@tools$IREA$table$effect)))
})

test_that("human orthologue mapping is explicit and conflicts are omitted", {
  r <- make_irea_fixture()
  r$species <- "Human"
  r$cytokine$Gene_Human <- c("H1","H4")
  r$polarization$Gene_Human <- "H4"
  out <- RunIREA(r,genes="H1")
  expect_equal(out$matched_genes,"G1")
  expect_equal(out$parameters$species,"Human")
  expect_error(RunIREA(r,genes="BAD"),"No unambiguous")
})

test_that("portal-style matrix file selects a named contrast", {
  r <- make_irea_fixture()
  path <- tempfile(fileext=".csv")
  on.exit(unlink(path),add=TRUE)
  utils::write.csv(data.frame(gene=c("G1","G2"),case1=c(1,-1),
                              case2=c(-1,1)),path,row.names=FALSE)
  expect_error(RunIREA(r,matrix=path,gene_diff_cutoff=0),"Select a contrast")
  a <- RunIREA(r,matrix=path,contrast="case1",gene_diff_cutoff=0)
  b <- RunIREA(r,matrix=path,contrast="case2",gene_diff_cutoff=0)
  expect_equal(a$table$effect,-b$table$effect)
})

test_that("no positively enriched polarization has a zero radar score", {
  r <- make_irea_fixture()
  out <- RunIREA(r,genes="G3",analysis="cell_polarization")
  expect_true(all(out$table$radar_score == 0))
  expect_true(all(out$table$effect <= 0))
  expect_error(PrepareIREAReference(tempdir(),"unknown"),"Unsupported")
})

test_that("multiple results can be compared without dropping terms", {
  r <- make_irea_fixture()
  a <- RunIREA(r,genes="G1")
  b <- RunIREA(r,genes="G4")
  plot <- IREAPlot(list(first=a,second=b),"heatmap")
  expect_s3_class(plot,"ggplot")
  expect_equal(nrow(plot$data),nrow(a$table)+nrow(b$table))
})

test_that("sparse and dense reference layers give the same scores", {
  sparse <- make_irea_fixture()
  dense <- sparse
  data <- as.matrix(SeuratObject::LayerData(dense$object,assay="RNA",layer="data"))
  SeuratObject::LayerData(dense$object,assay="RNA",layer="data") <- data
  a <- RunIREA(sparse,matrix=c(G1=1,G2=-1),gene_diff_cutoff=0)
  b <- RunIREA(dense,matrix=c(G1=1,G2=-1),gene_diff_cutoff=0)
  expect_equal(a$table$effect,b$table$effect)
  expect_equal(a$table$p_value,b$table$p_value)
})

test_that("polarization compares against PBS rather than unpolarized treated cells", {
  r <- make_irea_fixture()
  d <- SeuratObject::LayerData(r$object, assay="RNA", layer="data")
  score <- as.numeric(d["G1",])
  actual <- RunIREA(r,genes="G1",analysis="cell_polarization")$table
  target <- score[r$object$polarization == "NK-e"]
  control <- score[r$object$sample == "PBS"]
  expect_equal(actual$effect, mean(target)-mean(control))
  expect_equal(actual$p_value, stats::wilcox.test(target,control,exact=FALSE)$p.value)
  expect_equal(actual$n_control, 4L)
  expect_false(isTRUE(all.equal(actual$effect,
    mean(target)-mean(score[r$object$polarization == "None"]))))
  r$object$sample[r$object$sample == "PBS"] <- "Other"
  expect_error(RunIREA(r,genes="G1",analysis="cell_polarization"),"PBS baseline")
})

test_that("Seurat group labels are aligned to cells in the requested layer", {
  r <- make_irea_fixture()
  o <- r$object
  counts <- SeuratObject::LayerData(o, assay="RNA", layer="counts")
  selected <- c(1,2,5,6)
  o[["partial"]] <- SeuratObject::CreateAssay5Object(counts=counts[,selected])
  actual <- .irea_input(NULL,NULL,o,"sample","IL12","PBS","partial","counts",NULL)
  expected <- Matrix::rowMeans(counts[,5:6])-Matrix::rowMeans(counts[,1:2])
  expect_equal(actual$matrix, expected)
  expect_error(.irea_input(NULL,NULL,o,"sample","IL15","PBS","partial","counts",NULL),
               "selected assay/layer")
  o$sample[3] <- NA_character_
  actual <- .irea_input(NULL,NULL,o,"sample","IL12","PBS","RNA","counts",NULL)
  expect_equal(actual$matrix, Matrix::rowMeans(counts[,5:8])-Matrix::rowMeans(counts[,c(1,2,4)]))
  expect_error(RunIREA(r,object=o,group_by="sample",case="PBS",control="PBS"),"distinct")
})

test_that("comparison plots reject incompatible effect scales", {
  r <- make_irea_fixture()
  a <- RunIREA(r,genes="G1")
  b <- RunIREA(r,genes="G1",method="hypergeometric")
  expect_error(IREAPlot(list(a=a,b=b),"heatmap"),"same method")
  expect_error(IREAPlot(list(a=a,b=b),"dotplot"),"same method")
})

test_that("radar thresholds are recomputed without changing result statistics", {
  r <- RunIREA(make_irea_fixture(),genes="G4",analysis="cell_polarization")
  r$table <- data.frame(term=c("A","B","C"),effect=c(2,1,-1),
                        fdr=c(0.03,0.2,0.8),radar_score=c(1,0.5,0))
  before <- r$table
  loose <- IREAPlot(r,"radar",fdr_cutoff=0.05)$layers[[3]]$data
  strict <- IREAPlot(r,"radar",fdr_cutoff=0.01)$layers[[3]]$data
  expect_equal(loose$radar_score,c(1,0.5,0))
  expect_equal(strict$radar_score,c(0,0,0))
  expect_equal(r$table,before)
  expect_equal(.irea_radar_score(c(NA,1),c(NA,0.03)),c(0,1))
})

test_that("missing list genes and nonfinite cutoffs are handled explicitly", {
  r <- make_irea_fixture()
  expect_equal(RunIREA(r,genes=c(NA,"G1",""))$matched_genes,"G1")
  expect_error(RunIREA(r,genes=NA_character_),"No usable")
  expect_error(RunIREA(r,genes="G1",gene_diff_cutoff=Inf),"nonnegative")
})

test_that("split layer names cannot silently select the first batch", {
  r <- make_irea_fixture()
  counts <- SeuratObject::LayerData(r$object,assay="RNA",layer="counts")
  o <- SeuratObject::CreateSeuratObject(counts=list(first=counts[,1:6],second=counts[,7:12]))
  o$sample <- r$object$sample
  expect_error(.irea_input(NULL,NULL,o,"sample","IL12","PBS","RNA","counts",NULL),
               "exact name")
})

test_that("multiple matrix columns share the specified BH family", {
  r <- make_irea_fixture()
  input <- data.frame(gene=c("G1","G2","G3"),a=c(1,-1,0),b=c(0,0,1))
  a <- RunIREA(r,matrix=input,contrast="a",gene_diff_cutoff=0)
  b <- RunIREA(r,matrix=input,contrast="b",gene_diff_cutoff=0)
  expected <- stats::p.adjust(c(a$table$p_value,b$table$p_value),"BH")
  expect_equal(c(a$table$fdr,b$table$fdr),expected)
  expect_equal(a$parameters$contrasts_adjusted,c("a","b"))
  single <- RunIREA(r,matrix=input,contrast="a",gene_diff_cutoff=0,
                    fdr_scope="selected_contrast")
  expect_equal(single$table$fdr,stats::p.adjust(single$table$p_value,"BH"))
  expect_equal(single$table$p_value,a$table$p_value)
  expect_equal(single$table$effect,a$table$effect)
  expect_equal(a$parameters$fdr_family,"reference_terms_across_all_supplied_contrasts")
})

test_that("published NK polarization score statistics remain concordant", {
  fixture <- readRDS(test_path("fixtures","irea","nk_polarization_score.rds"))
  groups <- list(group=fixture$polarization,
                 terms=sort(setdiff(unique(fixture$polarization),"None")),
                 baseline=which(fixture$sample == "PBS"))
  actual <- .irea_score_table(fixture$score,groups)
  expected <- fixture$expected[match(actual$term,fixture$expected$Polarization),]
  expect_equal(actual$effect,expected[["Enrichment Score"]],tolerance=1e-8)
  expect_true(all(abs(actual$p_value/expected$pval-1) < 1e-8))
  expect_true(all(abs(actual$fdr/expected$padj-1) < 1e-8))
  # The earlier None comparator must not pass this external numerical fixture.
  groups$baseline <- which(fixture$polarization == "None")
  wrong <- .irea_score_table(fixture$score,groups)
  expect_gt(max(abs(wrong$effect-expected[["Enrichment Score"]])),0.4)
})
