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
