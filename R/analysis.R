searchDatasets <- function(datasets, pattern) {
  colIdx <- vapply(
    datasets,
    FUN = function(x) {
      return(
        grepl(
          pattern = pattern,
          x = x,
          ignore.case = TRUE
        )
      )
    },
    FUN.VALUE = logical(length = nrow(datasets))
  )
  rowIdx <- apply(colIdx, 1, any)
  if (any(rowIdx)) {
    return(datasets[rowIdx, , drop = FALSE])
  } else {
    message("No matching datasets found")
    return(NULL)
  }
}

#' Prefetch cycle gene
#'
#' Based on the human cell cycle genes, the cell cycle genes of the corresponding species were captured by homologous gene conversion.
#'
#' @inheritParams GeneConvert
#' @param species Latin names for animals,i.e., "Homo_sapiens", "Mus_musculus"
#' @param use_cached_gene Whether to use previously cached cell cycle gene conversion results for the species.
#'
#' @return A list of S-phase and G2M-phase genes.
#'
#' \code{\link{GeneConvert}}
#'
#' @export
#'
#' @examples
#' ccgenes <- CC_GenePrefetch("Homo_sapiens")
#' str(ccgenes)
#' ccgenes <- CC_GenePrefetch("Mus_musculus")
#' str(ccgenes)
CC_GenePrefetch <- function(
    species = "Homo_sapiens",
    Ensembl_version = 103,
    mirror = NULL,
    max_tries = 5,
    use_cached_gene = TRUE) {
  S <- Seurat::cc.genes.updated.2019$s.genes
  G2M <- Seurat::cc.genes.updated.2019$g2m.genes
  res <- NULL
  if (species != "Homo_sapiens") {
    if (isTRUE(use_cached_gene)) {
      res <- R.cache::loadCache(key = list(species))
    }
    if (is.null(res)) {
      res <- GeneConvert(
        geneID = unique(c(S, G2M)),
        geneID_from_IDtype = "symbol",
        geneID_to_IDtype = "symbol",
        species_from = "Homo_sapiens",
        species_to = species,
        Ensembl_version = Ensembl_version,
        max_tries = max_tries,
        mirror = mirror
      )
      R.cache::saveCache(res, key = list(species))
    } else {
      message("Using cached conversion results for ", species)
    }
    genes <- res[["geneID_collapse"]]
    S <- unlist(genes[S[S %in% rownames(genes)], "symbol"])
    G2M <- unlist(genes[G2M[G2M %in% rownames(genes)], "symbol"])
  }
  return(
    list(
      res = res,
      S = S,
      G2M = G2M
    )
  )
}

metap <- function(
    p,
    method = c(
      "maximump",
      "minimump",
      "wilkinsonp",
      "meanp",
      "sump",
      "votep"
    ),
    ...) {
  method <- match.arg(method)
  res <- do.call(method, args = list(p = p, ...))
  return(res)
}

wilkinsonp <- function(p, r = 1, alpha = 0.05, log.p = FALSE) {
  alpha <- ifelse(alpha > 1, alpha / 100, alpha)
  stopifnot(alpha > 0, alpha < 1)
  alpha <- ifelse(alpha > 0.5, 1 - alpha, alpha)
  keep <- (p >= 0) & (p <= 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    warning("Must have at least two valid p values")
    res <- list(
      p = NA_real_,
      pr = NA_real_,
      r = r,
      critp = NA_real_,
      alpha = alpha,
      validp = p[keep]
    )
  } else {
    pi <- p[keep]
    k <- length(pi)
    if (k != length(p)) {
      warning("Some studies omitted")
    }
    if ((r < 1) | (r > k)) {
      r <- 1
      warning("Illegal r set to 1")
    }
    pi <- sort(pi)
    pr <- pi[r]
    res <- list(
      p = stats::pbeta(pr, r, k + 1 - r, log.p = log.p),
      pr = pr,
      r = r,
      critp = stats::qbeta(alpha, r, k + 1 - r),
      alpha = alpha,
      validp = pi
    )
  }
  res
}

maximump <- function(p, alpha = 0.05, log.p = FALSE) {
  keep <- (p >= 0) & (p <= 1)
  validp <- p[keep]
  k <- length(validp)
  res <- wilkinsonp(p, r = k, alpha, log.p)
  res
}

minimump <- function(p, alpha = 0.05, log.p = FALSE) {
  res <- wilkinsonp(p, r = 1, alpha, log.p)
  res
}

meanp <- function(p) {
  keep <- (p >= 0) & (p <= 1)
  invalid <- sum(1L * keep) < 4
  if (invalid) {
    warning("Must have at least four valid p values")
    res <- list(z = NA_real_, p = NA_real_, validp = p[keep])
  } else {
    pi <- mean(p[keep])
    k <- length(p[keep])
    z <- (0.5 - pi) * sqrt(12 * k)
    if (k != length(p)) {
      warning("Some studies omitted")
    }
    res <- list(
      z = z,
      p = stats::pnorm(z, lower.tail = FALSE),
      validp = p[keep]
    )
  }
  res
}

sump <- function(p) {
  keep <- (p >= 0) & (p <= 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    warning("Must have at least two valid p values")
    res <- list(p = NA_real_, conservativep = NA_real_, validp = p[keep])
  } else {
    sigmap <- sum(p[keep])
    k <- length(p[keep])
    conservativep <- exp(k * log(sigmap) - lgamma(k + 1))
    nterm <- floor(sigmap) + 1
    denom <- lfactorial(k)
    psum <- 0
    terms <- vector("numeric", nterm)
    for (i in 1:nterm) {
      terms[i] <- lchoose(k, i - 1) +
        k *
          log(
            sigmap -
              i +
              1
          ) -
        denom
      pm <- 2 * (i %% 2) - 1
      psum <- psum + pm * exp(terms[i])
    }
    if (k != length(p)) {
      warning("Some studies omitted")
    }
    if (sigmap > 20) {
      warning("Likely to be unreliable, check with another method")
    }
    res <- list(
      p = psum,
      conservativep = conservativep,
      validp = p[keep]
    )
  }
  res
}

votep <- function(p, alpha = 0.5) {
  alpha <- ifelse(alpha > 1, alpha / 100, alpha)
  stopifnot(alpha > 0, alpha < 1)
  keep <- (p >= 0) & (p <= 1)
  alp <- vector("numeric", 2)
  if (alpha <= 0.5) {
    alp[1] <- alpha
    alp[2] <- 1 - alpha
  } else {
    alp[2] <- alpha
    alp[1] <- 1 - alpha
  }
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    warning("Must have at least two valid p values")
    res <- list(
      p = NA_real_,
      pos = NA_integer_,
      neg = NA_integer_,
      alpha = alpha,
      validp = p[keep]
    )
  } else {
    pi <- p[keep]
    k <- length(pi)
    pos <- sum(1L * (pi < alp[1]))
    neg <- sum(1L * (pi > alp[2]))
    if (k != length(p)) {
      warning("Some studies omitted")
    }
    if ((pos + neg) <= 0) {
      warning("All p values are within specified limits of alpha")
      p <- 1
    } else {
      p <- stats::binom.test(
        pos, pos + neg, 0.5,
        alternative = "greater"
      )$p.value
    }
    res <- list(
      p = p,
      pos = pos,
      neg = neg,
      alpha = alpha,
      validp = pi
    )
  }
  res
}

py_to_r_auto <- function(x) {
  if (inherits(x, "python.builtin.object")) {
    x <- reticulate::py_to_r(x)
  }
  return(x)
}

#' Convert a seurat object to an anndata object using reticulate
#'
#' This function takes a Seurat object and converts it to an anndata object using the reticulate package.
#'
#' @param srt A Seurat object.
#' @param assay_x Assay to convert as the main data matrix (X) in the anndata object.
#' @param layer_x Layer name for assay_x in the Seurat object.
#' @param assay_y Assays to convert as layers in the anndata object.
#' @param layer_y Layer names for the assay_y in the Seurat object.
#' @param convert_tools Logical indicating whether to convert the tool-specific data.
#' @param convert_misc Logical indicating whether to convert the miscellaneous data.
#' @param features Optional vector of features to include in the anndata object.
#' Defaults to all features in assay_x.
#' @param verbose Logical indicating whether to print verbose messages during the conversion process.
#'
#' @return A \code{anndata} object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
#' adata <- srt_to_adata(pancreas_sub)
#' adata
#'
#' ## Or save as an h5ad file or a loom file
#' # adata$write_h5ad(
#' #   "pancreas_sub.h5ad"
#' # )
#' # adata$write_loom(
#' #   "pancreas_sub.loom",
#' #   write_obsm_varm = TRUE
#' # )
#' }
srt_to_adata <- function(
    srt,
    features = NULL,
    assay_x = "RNA",
    layer_x = "counts",
    assay_y = c("spliced", "unspliced"),
    layer_y = "counts",
    convert_tools = FALSE,
    convert_misc = FALSE,
    verbose = TRUE) {
  check_python(c("scanpy", "numpy"))

  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }

  if (is.null(features)) {
    features <- rownames(srt[[assay_x]])
  }
  if (length(layer_y) == 1) {
    layer_y <- rep(layer_y, length(assay_y))
    names(layer_y) <- assay_y
  } else if (length(layer_y) != length(assay_y)) {
    stop("layer_y must be one character or the same length of the assay_y")
  }

  sc <- reticulate::import("scanpy", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  obs <- srt@meta.data
  if (ncol(obs) > 0) {
    for (i in seq_len(ncol(obs))) {
      if (is.logical(obs[, i])) {
        obs[, i] <- factor(
          as.character(obs[, i]),
          levels = c("TRUE", "FALSE")
        )
      }
    }
  }

  var <- GetFeaturesData(srt, assay = assay_x)[features, , drop = FALSE]
  if (ncol(var) > 0) {
    for (i in seq_len(ncol(var))) {
      if (
        is.logical(var[, i]) && !identical(colnames(var)[i], "highly_variable")
      ) {
        var[, i] <- factor(
          as.character(var[, i]),
          levels = c("TRUE", "FALSE")
        )
      }
    }
  }
  if (length(SeuratObject::VariableFeatures(srt, assay = assay_x) > 0)) {
    if ("highly_variable" %in% colnames(var)) {
      var <- var[, colnames(var) != "highly_variable"]
    }
    var[["highly_variable"]] <- features %in%
      SeuratObject::VariableFeatures(srt, assay = assay_x)
  }

  X <- Matrix::t(
    SeuratObject::GetAssayData(
      srt,
      assay = assay_x,
      layer = layer_x
    )[features, , drop = FALSE]
  )
  adata <- sc$AnnData(
    X = reticulate::np_array(X, dtype = np$float32),
    obs = obs,
    var = cbind(
      data.frame(features = features),
      var
    )
  )
  adata$var_names <- features

  layer_list <- list()
  for (assay in names(srt@assays)[names(srt@assays) != assay_x]) {
    if (assay %in% assay_y) {
      layer <- Matrix::t(
        SeuratObject::GetAssayData(
          srt,
          assay = assay,
          layer = layer_y[assay]
        )
      )
      if (!identical(dim(layer), dim(X))) {
        if (all(colnames(X) %in% colnames(layer))) {
          layer <- layer[, colnames(X)]
        } else {
          stop(
            "The following features in the '",
            assay_x,
            "' assay can not be found in the '",
            assay,
            "' assay:\n  ",
            paste0(
              utils::head(colnames(X)[!colnames(X) %in% colnames(layer)], 10),
              collapse = ","
            ),
            "..."
          )
        }
      }
      layer_list[[assay]] <- layer
    } else {
      if (isTRUE(verbose)) {
        message("Assay '", assay, "' is in the srt object but not converted.")
      }
    }
  }
  if (length(layer_list) > 0) {
    adata$layers <- layer_list
  }

  reduction_list <- list()
  for (reduction in names(srt@reductions)) {
    reduction_list[[paste0(reduction)]] <- srt[[reduction]]@cell.embeddings
  }
  if (length(reduction_list) > 0) {
    adata$obsm <- reduction_list
  }

  obsp_list <- list()
  for (graph in names(srt@graphs)) {
    obsp_list[[graph]] <- srt[[graph]]
  }
  for (neighbor in names(srt@neighbors)) {
    obsp_list[[neighbor]] <- srt[[neighbor]]
  }
  if (length(obsp_list) > 0) {
    adata$obsp <- obsp_list
  }

  uns_list <- list()
  if (isTRUE(convert_misc)) {
    for (nm in names(srt@misc)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@misc[[nm]]
      }
    }
  } else {
    if (isTRUE(verbose)) {
      message("'misc' slot is not converted.")
    }
  }
  if (isTRUE(convert_tools)) {
    for (nm in names(srt@tools)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@tools[[nm]]
      }
    }
  } else {
    if (isTRUE(verbose)) {
      message("'tools' slot is not converted.")
    }
  }
  if (length(uns_list) > 0) {
    adata$uns <- uns_list
  }

  return(adata)
}

max_depth <- function(x, depth = 0) {
  if (is.list(x)) {
    return(max(unlist(lapply(x, max_depth, depth + 1))))
  } else {
    return(depth)
  }
}

check_python_element <- function(x, depth = max_depth(x)) {
  if (depth == 0 || !is.list(x) || !inherits(x, "python.builtin.object")) {
    if (inherits(x, "python.builtin.object")) {
      x_r <- tryCatch(reticulate::py_to_r(x), error = identity)
      if (inherits(x_r, "error")) {
        return(x)
      } else {
        return(x_r)
      }
    } else {
      return(x)
    }
  } else {
    raw_depth <- max_depth(x)
    x <- lapply(x, function(element) {
      if (inherits(element, "python.builtin.object")) {
        element_r <- tryCatch(reticulate::py_to_r(element), error = identity)
        if (inherits(element_r, "error")) {
          return(element)
        } else {
          return(element_r)
        }
      } else {
        return(element)
      }
    })
    cur_depth <- max_depth(x)
    if (cur_depth > raw_depth) {
      depth <- depth + 1
    }
    x_checked <- lapply(x, check_python_element, depth - 1)
    return(x_checked)
  }
}
