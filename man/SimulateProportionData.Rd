% Generated manually because roxygen2 is unavailable in this environment
% Please edit documentation in R/SimulateProportionData.R
\name{SimulateProportionData}
\alias{SimulateProportionData}
\title{Simulate Single-Cell Data for Proportion Testing}
\usage{
SimulateProportionData(
  n_genes = 2000L,
  samples_per_condition = 5L,
  cells_per_sample = 350L,
  conditions = c("Control", "Treatment"),
  cell_types = c(
    "Tcell",
    "Bcell",
    "Myeloid",
    "NK",
    "DC",
    "Mono",
    "Epithelial",
    "Fibro",
    "Mast"
  ),
  base_props = c(0.2, 0.14, 0.16, 0.1, 0.08, 0.1, 0.12, 0.07, 0.03),
  condition_prop_shift = c(0.07, -0.03, -0.03, -0.01, -0.01, -0.01, 0.01, 0,
    0.01),
  condition_prop_matrix = NULL,
  prop_concentration = 35,
  marker_genes_per_type = 30L,
  marker_fc = 4,
  condition_de_genes = 100L,
  condition_de_fc = 1.5,
  dispersion = 2,
  sample_prefix = "S",
  batch_levels = c("Batch1", "Batch2", "Batch3"),
  assay = "RNA",
  seed = 11,
  verbose = TRUE
)
}
\arguments{
\item{n_genes}{Number of genes to simulate.}
\item{samples_per_condition}{Number of biological samples per condition.}
\item{cells_per_sample}{Number of cells per sample.}
\item{conditions}{Condition labels.}
\item{cell_types}{Cell type labels.}
\item{base_props}{Baseline cell-type proportions for the first condition.}
\item{condition_prop_shift}{Additive shift applied to baseline proportions for the second condition.}
\item{condition_prop_matrix}{Optional condition-specific proportion matrix (rows=conditions, cols=cell types).}
\item{prop_concentration}{Dirichlet concentration controlling sample-to-sample variability.}
\item{marker_genes_per_type}{Number of marker genes for each cell type.}
\item{marker_fc}{Fold-change multiplier for marker genes.}
\item{condition_de_genes}{Number of global condition-affected genes.}
\item{condition_de_fc}{Fold-change multiplier for global condition-affected genes.}
\item{dispersion}{Negative binomial size parameter.}
\item{sample_prefix}{Prefix for sample IDs.}
\item{batch_levels}{Batch labels assigned to samples.}
\item{assay}{Assay name for Seurat object creation.}
\item{seed}{Random seed.}
\item{verbose}{Whether to print messages.}
}
\description{
Generate a synthetic \code{Seurat} object with sample-level replicates and
condition-driven composition shifts, designed for DA testing workflows.
Default settings generate richer cell-type diversity (including rare groups)
and sample structure for more realistic proportion-test regression checks.
}
\value{
A \code{Seurat} object containing metadata columns:
\code{Sample}, \code{Condition}, \code{Batch}, \code{CellType}, \code{SubCellType}, and \code{Phase}.
Simulation truth is stored in \code{srt@tools[["SimulateProportionData"]]}.
}
\examples{
\dontrun{
srt_sim <- SimulateProportionData(
  n_genes = 1500,
  samples_per_condition = 4,
  cells_per_sample = 250,
  seed = 42
)

srt_sim <- RunProportionTest(
  srt_sim,
  group.by = "CellType",
  split.by = "Condition",
  sample.by = "Sample",
  proportion_method = "propeller",
  comparison = list(c("Control", "Treatment"))
)
}
}
\seealso{\link{RunProportionTest}}
