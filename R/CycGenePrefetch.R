#' @title Prefetch cell cycle genes
#'
#' @description
#' Based on the human cell cycle genes,
#' the cell cycle genes of the corresponding species were captured by homologous gene conversion.
#'
#' @md
#' @inheritParams GeneConvert
#' @inheritParams thisutils::log_message
#' @param species Latin names for animals, i.e., `"Homo_sapiens"`, `"Mus_musculus"`
#' @param use_cached_gene Whether to use previously cached cell cycle gene conversion results for the species.
#' Default is `TRUE`.
#'
#' @return A list of S-phase and G2M-phase genes.
#'
#' @seealso [GeneConvert]
#'
#' @export
#'
#' @examples
#' ccgenes <- CycGenePrefetch("Homo_sapiens")
#' str(ccgenes)
#'
#' ccgenes <- CycGenePrefetch("Mus_musculus")
#' str(ccgenes)
CycGenePrefetch <- function(
    species = "Homo_sapiens",
    Ensembl_version = NULL,
    mirror = NULL,
    max_tries = 5,
    use_cached_gene = TRUE,
    verbose = TRUE) {
  log_message(
    "Prefetching cell cycle genes for {.val {species}} ...",
    verbose = verbose
  )
  s_genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes
  res <- NULL
  if (species != "Homo_sapiens") {
    if (isTRUE(use_cached_gene)) {
      check_r("R.cache", verbose = FALSE)
      res <- R.cache::loadCache(key = list(species))
    }
    if (is.null(res)) {
      res <- GeneConvert(
        geneID = unique(c(s_genes, g2m_genes)),
        geneID_from_IDtype = "symbol",
        geneID_to_IDtype = "symbol",
        species_from = "Homo_sapiens",
        species_to = species,
        Ensembl_version = Ensembl_version,
        max_tries = max_tries,
        mirror = mirror,
        verbose = verbose
      )
      R.cache::saveCache(res, key = list(species))
      log_message(
        "Cached conversion results for {.val {species}}",
        verbose = verbose
      )
    } else {
      log_message(
        "Using cached conversion results for {.val {species}}",
        verbose = verbose
      )
    }
    genes <- res[["geneID_collapse"]]
    s_genes <- unlist(
      genes[s_genes[s_genes %in% rownames(genes)], "symbol"]
    )
    g2m_genes <- unlist(
      genes[g2m_genes[g2m_genes %in% rownames(genes)], "symbol"]
    )
  }
  log_message(
    "Cell cycle gene prefetching completed {.val {species}}",
    message_type = "success",
    verbose = verbose
  )
  return(
    list(
      res = res,
      S = s_genes,
      G2M = g2m_genes
    )
  )
}
