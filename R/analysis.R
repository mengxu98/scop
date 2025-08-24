#' @title Prefetch cell cycle genes
#'
#' @description
#' Based on the human cell cycle genes, the cell cycle genes of the corresponding species were captured by homologous gene conversion.
#'
#' @md
#' @inheritParams GeneConvert
#' @param species Latin names for animals,i.e., `"Homo_sapiens"`, `"Mus_musculus"`
#' @param use_cached_gene Whether to use previously cached cell cycle gene conversion results for the species.
#' @param verbose Whether to print messages.
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
    "Prefetching cell cycle genes for {.val {species}}...",
    verbose = verbose
  )
  s_genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes
  res <- NULL
  if (species != "Homo_sapiens") {
    if (isTRUE(use_cached_gene)) {
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
        mirror = mirror
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

wilkinsonp <- function(p, r = 1, alpha = 0.05, log.p = FALSE) {
  alpha <- ifelse(alpha > 1, alpha / 100, alpha)
  stopifnot(alpha > 0, alpha < 1)
  alpha <- ifelse(alpha > 0.5, 1 - alpha, alpha)
  keep <- (p >= 0) & (p <= 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    log_message(
      "Must have at least two valid p values",
      message_type = "warning"
    )
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
      log_message(
        "Some studies omitted",
        message_type = "warning"
      )
    }
    if ((r < 1) | (r > k)) {
      r <- 1
      log_message(
        "Illegal r set to 1",
        message_type = "warning"
      )
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
    log_message(
      "Must have at least four valid p values",
      message_type = "warning"
    )
    res <- list(z = NA_real_, p = NA_real_, validp = p[keep])
  } else {
    pi <- mean(p[keep])
    k <- length(p[keep])
    z <- (0.5 - pi) * sqrt(12 * k)
    if (k != length(p)) {
      log_message(
        "Some studies omitted",
        message_type = "warning"
      )
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
    log_message(
      "Must have at least two valid p values",
      message_type = "warning"
    )
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
      log_message(
        "Some studies omitted",
        message_type = "warning"
      )
    }
    if (sigmap > 20) {
      log_message(
        "Likely to be unreliable, check with another method",
        message_type = "warning"
      )
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
    log_message(
      "Must have at least two valid p values",
      message_type = "warning"
    )
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
      log_message(
        "Some studies omitted",
        message_type = "warning"
      )
    }
    if ((pos + neg) <= 0) {
      log_message(
        "All p-values are within specified limits of alpha",
        message_type = "warning"
      )
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

max_depth <- function(x, depth = 0) {
  if (is.list(x)) {
    return(max(unlist(lapply(x, max_depth, depth + 1))))
  } else {
    return(depth)
  }
}

check_python_element <- function(
    x,
    depth = max_depth(x)) {
  if (depth == 0 || !is.list(x) || !inherits(x, "python.builtin.object")) {
    if (inherits(x, "python.builtin.object")) {
      x_r <- tryCatch(py_to_r2(x), error = identity)
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
        element_r <- tryCatch(py_to_r2(element), error = identity)
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
