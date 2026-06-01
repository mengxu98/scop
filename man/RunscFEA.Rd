\name{RunscFEA}
\alias{RunscFEA}
\title{Run scFEA flux estimation for a Seurat object}
\usage{
RunscFEA(
  srt,
  assay = NULL,
  layer = "data",
  species = c("human", "mouse"),
  n_epoch = 100,
  sc_imputation = FALSE,
  assay_flux = "scFEAflux",
  assay_balance = "scFEAbalance",
  store_metadata = FALSE,
  data_dir = NULL,
  seed = 16,
  verbose = TRUE
)
}
\arguments{
\item{srt}{A Seurat object.}
\item{assay}{Assay to use as expression matrix. Default is \code{DefaultAssay(srt)}.}
\item{layer}{Assay layer to use. Default is \code{"data"}.}
\item{species}{One of \code{"human"} or \code{"mouse"}, selecting the M168 scFEA files.}
\item{n_epoch}{Number of scFEA training epochs.}
\item{sc_imputation}{Whether to run MAGIC imputation inside the scFEA backend.}
\item{assay_flux}{Name of the assay storing module flux scores.}
\item{assay_balance}{Name of the assay storing metabolite balance scores.}
\item{store_metadata}{Whether to also append flux and balance values to \code{srt@meta.data}.}
\item{data_dir}{Optional directory containing scFEA M168 CSV resources. If \code{NULL}, files are downloaded from \code{mengxu98/datasets} and cached with \code{tools::R_user_dir("scop", "data")}.}
\item{seed}{Random seed passed to R and the Python scFEA backend.}
\item{verbose}{Whether to print progress messages.}
}
\value{
A Seurat object with \code{assay_flux}, \code{assay_balance}, and \code{srt@tools[["scFEA"]]}.
}
\description{
\code{RunscFEA()} calls the bundled Python scFEA backend through \pkg{reticulate}. The scFEA GNN architecture, loss function, and training loop are kept in the backend. M168 human/mouse model files are downloaded from the \code{mengxu98/datasets} repository and cached outside the package; the R side prepares a Seurat expression matrix and stores flux / balance outputs back into Seurat.
}
\details{
scFEA is licensed for academic, non-commercial use. See the bundled \code{inst/python/scfea/LICENSE} file for details. Data resources are downloaded from \url{https://github.com/mengxu98/datasets/tree/main/scFEA}.
}
