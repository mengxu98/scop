#' Convert homologous gene symbols in expression objects
#'
#' @description
#' Convert feature names between species with [GeneConvert] and collapse
#' duplicated target homologs by summing expression values. The Seurat method
#' rebuilds the selected assay from the converted counts matrix and keeps cell
#' metadata and spatial images when present.
#'
#' @md
#' @inheritParams GeneConvert
#' @param object A `Seurat` object or a gene-by-cell matrix.
#' @param assay Assay to convert when `object` is a `Seurat` object. If `NULL`,
#' the default assay is used.
#' @param layer Assay layer used for conversion. Default `"counts"`.
#' @param multi_mapping How to handle source genes mapped to multiple target
#' homologs. `"first"` keeps the first target homolog for each source gene.
#' @param keep_unmapped Whether to keep unmapped source genes with their
#' original names.
#' @param collapse_fun Function used to collapse duplicated target homologs.
#' Currently only `"sum"` is supported.
#'
#' @return A converted object of the same high-level type as `object`. The
#' mapping table is stored in `@tools$ConvertHomologs` for Seurat objects and in
#' the `"ConvertHomologs"` attribute for matrix inputs.
#'
#' @seealso
#' [AnnotateFeatures], [ConvertHomologs]]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_human <- ConvertHomologs(
#'   pancreas_sub,
#'   species_from = "Mus_musculus",
#'   species_to = "Homo_sapiens"
#' )
#' rownames(pancreas_human)[1:5]
ConvertHomologs <- function(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
) {
  check_r("biomaRt", verbose = FALSE)
  if (inherits(object, "Seurat")) {
    return(ConvertHomologs.Seurat(
      object = object,
      species_from = species_from,
      species_to = species_to,
      geneID_from_IDtype = geneID_from_IDtype,
      geneID_to_IDtype = geneID_to_IDtype,
      assay = assay,
      layer = layer,
      multi_mapping = multi_mapping,
      keep_unmapped = keep_unmapped,
      collapse_fun = collapse_fun,
      Ensembl_version = Ensembl_version,
      biomart = biomart,
      mirror = mirror,
      max_tries = max_tries,
      verbose = verbose
    ))
  }
  if (is.matrix(object) || inherits(object, "Matrix")) {
    return(ConvertHomologs.matrix(
      object = object,
      species_from = species_from,
      species_to = species_to,
      geneID_from_IDtype = geneID_from_IDtype,
      geneID_to_IDtype = geneID_to_IDtype,
      assay = assay,
      layer = layer,
      multi_mapping = multi_mapping,
      keep_unmapped = keep_unmapped,
      collapse_fun = collapse_fun,
      Ensembl_version = Ensembl_version,
      biomart = biomart,
      mirror = mirror,
      max_tries = max_tries,
      verbose = verbose
    ))
  }
  UseMethod("ConvertHomologs")
}

#' @rdname ConvertHomologs
#' @export
ConvertHomologs.Seurat <- function(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  counts <- GetAssayData5(
    object = object,
    assay = assay,
    layer = layer
  )
  converted <- ConvertHomologs.default(
    object = counts,
    species_from = species_from,
    species_to = species_to,
    geneID_from_IDtype = geneID_from_IDtype,
    geneID_to_IDtype = geneID_to_IDtype,
    multi_mapping = multi_mapping,
    keep_unmapped = keep_unmapped,
    collapse_fun = collapse_fun,
    Ensembl_version = Ensembl_version,
    biomart = biomart,
    mirror = mirror,
    max_tries = max_tries,
    verbose = verbose
  )

  converted_srt <- Seurat::CreateSeuratObject(
    counts = converted,
    assay = assay,
    meta.data = object[[]],
    project = object@project.name
  )
  SeuratObject::DefaultAssay(converted_srt) <- assay
  converted_srt@images <- object@images
  converted_srt@misc <- object@misc
  converted_srt@tools <- object@tools
  converted_srt@tools$ConvertHomologs <- attr(converted, "ConvertHomologs")
  converted_srt
}

#' @rdname ConvertHomologs
#' @export
ConvertHomologs.matrix <- function(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
) {
  ConvertHomologs.default(
    object = object,
    species_from = species_from,
    species_to = species_to,
    geneID_from_IDtype = geneID_from_IDtype,
    geneID_to_IDtype = geneID_to_IDtype,
    multi_mapping = multi_mapping,
    keep_unmapped = keep_unmapped,
    collapse_fun = collapse_fun,
    Ensembl_version = Ensembl_version,
    biomart = biomart,
    mirror = mirror,
    max_tries = max_tries,
    verbose = verbose
  )
}

#' @rdname ConvertHomologs
#' @export
ConvertHomologs.Matrix <- ConvertHomologs.matrix

#' @rdname ConvertHomologs
#' @export
ConvertHomologs.default <- function(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
) {
  multi_mapping <- match.arg(multi_mapping)
  collapse_fun <- match.arg(collapse_fun)
  if (is.null(rownames(object))) {
    log_message(
      "{.arg object} must have feature names in {.code rownames(object)}",
      message_type = "error"
    )
  }

  conv <- GeneConvert(
    geneID = rownames(object),
    geneID_from_IDtype = geneID_from_IDtype,
    geneID_to_IDtype = geneID_to_IDtype,
    species_from = species_from,
    species_to = species_to,
    Ensembl_version = Ensembl_version,
    biomart = biomart,
    mirror = mirror,
    max_tries = max_tries,
    verbose = verbose
  )

  mapping <- ConvertHomologs_mapping(
    conv = conv,
    geneID_to_IDtype = geneID_to_IDtype,
    multi_mapping = multi_mapping,
    keep_unmapped = keep_unmapped,
    features = rownames(object)
  )
  if (nrow(mapping) == 0L) {
    log_message(
      "No homologous genes remained after conversion",
      message_type = "error"
    )
  }

  converted <- ConvertHomologs_collapse_matrix(
    mat = object[mapping$from_geneID, , drop = FALSE],
    target = mapping$to_geneID
  )
  attr(converted, "ConvertHomologs") <- list(
    mapping = mapping,
    species_from = species_from,
    species_to = species_to,
    geneID_from_IDtype = geneID_from_IDtype,
    geneID_to_IDtype = geneID_to_IDtype,
    Ensembl_version = conv$Ensembl_version,
    unmapped = conv$geneID_unmapped
  )
  log_message(
    "Converted {.val {length(unique(mapping$from_geneID))}} source genes to {.val {nrow(converted)}} target homologs",
    verbose = verbose
  )
  converted
}

ConvertHomologs_mapping <- function(
  conv,
  geneID_to_IDtype,
  multi_mapping,
  keep_unmapped,
  features
) {
  mapping <- data.frame(
    from_geneID = character(),
    to_geneID = character(),
    stringsAsFactors = FALSE
  )
  expand <- conv$geneID_expand
  if (!is.null(expand)) {
    expand <- as.data.frame(expand, stringsAsFactors = FALSE)
    to_col <- as.character(geneID_to_IDtype[[1]])
    if (!to_col %in% colnames(expand)) {
      fallback <- setdiff(
        colnames(expand),
        c("from_IDtype", "from_geneID", "to_IDtype", "to_geneID")
      )
      if (length(fallback) > 0L) {
        to_col <- fallback[[1]]
      }
    }
    if (all(c("from_geneID", to_col) %in% colnames(expand))) {
      mapping <- expand[, c("from_geneID", to_col), drop = FALSE]
      mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
      if (ncol(mapping) >= 2L) {
        mapping <- mapping[, seq_len(2L), drop = FALSE]
        colnames(mapping) <- c("from_geneID", "to_geneID")
        if (is.list(mapping$to_geneID)) {
          mapping <- do.call(
            rbind,
            Map(
              f = function(from, to) {
                to <- unique(as.character(to))
                to <- to[!is.na(to) & nzchar(to)]
                if (length(to) == 0L) {
                  return(NULL)
                }
                data.frame(
                  from_geneID = from,
                  to_geneID = to,
                  stringsAsFactors = FALSE
                )
              },
              mapping$from_geneID,
              mapping$to_geneID
            )
          )
          if (is.null(mapping)) {
            mapping <- data.frame(
              from_geneID = character(),
              to_geneID = character(),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  if (nrow(mapping) > 0L) {
    mapping$from_geneID <- as.character(mapping$from_geneID)
    mapping$to_geneID <- as.character(mapping$to_geneID)
    mapping <- mapping[
      !is.na(mapping$to_geneID) & nzchar(mapping$to_geneID), ,
      drop = FALSE
    ]
    mapping <- mapping[!duplicated(mapping), , drop = FALSE]
    if (identical(multi_mapping, "first")) {
      mapping <- mapping[!duplicated(mapping$from_geneID), , drop = FALSE]
    }
  }
  if (isTRUE(keep_unmapped)) {
    mapped <- unique(mapping$from_geneID)
    unmapped <- setdiff(features, mapped)
    if (length(unmapped) > 0L) {
      mapping <- rbind(
        mapping,
        data.frame(
          from_geneID = unmapped,
          to_geneID = unmapped,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  mapping
}

ConvertHomologs_collapse_matrix <- function(mat, target) {
  target <- as.character(target)
  target_levels <- unique(target)
  group_index <- match(target, target_levels)
  collapse_matrix <- Matrix::sparseMatrix(
    i = group_index,
    j = seq_along(target),
    x = 1,
    dims = c(length(target_levels), length(target)),
    dimnames = list(target_levels, rownames(mat))
  )
  converted <- collapse_matrix %*% mat
  methods::as(converted, "dgCMatrix")
}
