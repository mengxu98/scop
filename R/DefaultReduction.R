
#' @title Find the default reduction name in a Seurat object
#'
#' @param srt A Seurat object.
#' @param pattern Character string containing a regular expression to search for.
#' @param min_dim Minimum dimension threshold.
#' @param max_distance Maximum distance allowed for a match.
#'
#' @return Default reduction name.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' names(pancreas_sub@reductions)
#' DefaultReduction(pancreas_sub)
#'
#' # Searches for matches to "pca"
#' DefaultReduction(pancreas_sub, pattern = "pca")
#'
#' # Searches for approximate matches to "pc"
#' DefaultReduction(pancreas_sub, pattern = "pc")
DefaultReduction <- function(
    srt,
    pattern = NULL,
    min_dim = 2,
    max_distance = 0.1) {
  if (length(srt@reductions) == 0) {
    log_message(
      "Unable to find any reductions",
      message_type = "error"
    )
  }
  pattern_default <- c(
    "umap",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr",
    "pca",
    "svd",
    "ica",
    "nmf",
    "mds",
    "glmpca"
  )
  pattern_dim <- c("2D", "3D")
  reduc_all <- names(srt@reductions)
  reduc_all <- reduc_all[unlist(lapply(reduc_all, function(x) {
    dim(srt@reductions[[x]]@cell.embeddings)[2] >= min_dim
  }))]
  if (length(reduc_all) == 0) {
    log_message(
      "No dimensional reduction found in {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (length(reduc_all) == 1) {
    return(reduc_all)
  }
  if (is.null(pattern)) {
    if (("Default_reduction" %in% names(srt@misc))) {
      pattern <- srt@misc[["Default_reduction"]]
    } else {
      pattern <- pattern_default
    }
  }

  pattern <- c(pattern, paste0(pattern, min_dim, "D"))
  if (any(pattern %in% reduc_all)) {
    return(pattern[pattern %in% reduc_all][1])
  }
  index <- c(unlist(sapply(pattern, function(pat) {
    grep(pattern = pat, x = reduc_all, ignore.case = TRUE)
  })))
  if (length(index) > 0) {
    default_reduc <- reduc_all[index]
  } else {
    index <- c(
      unlist(
        sapply(
          pattern, function(pat) {
            agrep(
              pattern = pat,
              x = reduc_all,
              max.distance = max_distance,
              ignore.case = TRUE
            )
          }
        )
      )
    )
    if (length(index) > 0) {
      default_reduc <- reduc_all[index]
    } else {
      default_reduc <- reduc_all
    }
  }
  if (length(default_reduc) > 1) {
    default_reduc <- default_reduc[unlist(sapply(
      c(pattern_default, pattern_dim),
      function(pat) {
        grep(pattern = pat, x = default_reduc, ignore.case = TRUE)
      }
    ))]
    default_reduc <- default_reduc[which.min(sapply(
      default_reduc,
      function(x) dim(srt@reductions[[x]])[2]
    ))]
  }
  return(default_reduc)
}
