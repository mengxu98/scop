\name{GRNPlot}
\alias{GRNPlot}
\title{Plot networks from a unified GRN result}
\usage{
GRNPlot(
  object,
  group = NULL,
  layer = "tf_gene",
  type = c("network", "layered", "overlap"),
  top_n = 100,
  layout = "auto",
  label_nodes = TRUE,
  compare = NULL,
  compare_group = NULL,
  tool_name = NULL,
  title = NULL
)
}
\arguments{
\item{object}{A \code{GRNResult}, or a Seurat object containing one in \code{tools}.}
\item{group}{One result group. Required when the result contains multiple groups.}
\item{layer}{Network layer to display or compare.}
\item{type}{A network, a layered view of native chains, or a pairwise edge-overlap heatmap.}
\item{top_n}{Maximum edges in a network view, native chains in a layered view, or edges per method within the shared comparison scope. Use \code{Inf} for all entries.}
\item{layout}{One of \code{"auto"}, \code{"star"}, \code{"hub"}, \code{"tripartite"}, \code{"kk"} or \code{"fr"}. The automatic layout uses tripartite columns for layered views and Kamada--Kawai for other networks.}
\item{label_nodes}{Whether to label displayed nodes.}
\item{compare}{For \code{type = "overlap"}, another result or a named list of results. List names provide plot labels.}
\item{compare_group}{Group in each comparison result. Defaults to \code{group}. A vector supplies one group per comparison result.}
\item{tool_name}{Name of the Seurat \code{tools} entry. If omitted, the only
stored GRN result is used. Specify a name when several results are stored.
An overlap comparison uses the same name if present; otherwise it selects a
single stored result. Use \code{GetGRN(what = "result")} to select a differently
named result from a comparison with multiple networks.}
\item{title}{Optional plot title.}
}
\value{
A ggplot object. Its \code{grn_data} attribute contains the displayed edges and nodes, or overlap counts and the comparison scope. The original result is unchanged.
}
\description{
Read network views through \code{\link{GetGRN}} and draw nodes, edges and native chains while retaining direction, sign and membership evidence.
}
\details{
Signs are labelled positive, negative, zero or unknown. Nonnegative importance scores do not acquire an activation label. Membership edges have a separate evidence style. Arrowheads are drawn only for directed edges. Node shapes distinguish TFs, genes and regions.

A layered view selects existing native chains and then selects the corresponding TF--region and region--gene edges from their native layers. Each segment retains its own sign. Shared nodes can connect several selected chains; the displayed chain table is retained in the plot attribute.

Ranking and \code{top_n} apply only to the displayed view. Captions report the selected size. An overlap heatmap restricts each result to exactly matched source and target identifiers present in all compared results, then selects the top edges independently within each method. TF and gene endpoints share a gene identifier space. Directed and undirected edges remain distinct. Reversed endpoints of an undirected edge denote the same edge. Declared \code{input_summary$entity_namespace} values must agree. No identifier conversion, score aggregation or consensus inference is performed.

The network layout reuses the existing SCENIC graph-layout implementation. Layout and label placement are reproducible and preserve the caller's random-number state. Use existing \code{CellScoring}, \code{FeatureDimPlot} and \code{GroupHeatmap} functions to display extracted gene-set scores.
}
\seealso{\code{\link{RunGRN}}, \code{\link{GetGRN}}}
