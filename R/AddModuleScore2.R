#' @importFrom Seurat UpdateSymbolList CaseMatch
#' @importFrom SeuratObject DefaultAssay GetAssayData CheckGC
#' @importFrom BiocParallel bplapply bpaggregate
#' @importFrom stats rnorm
#' @importFrom Matrix rowMeans colMeans
#' @importFrom ggplot2 cut_number
AddModuleScore2 <- function(
    object,
    layer = "data",
    features,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    k = FALSE,
    assay = NULL,
    name = "Cluster",
    seed = 1,
    search = FALSE,
    BPPARAM = BiocParallel::bpparam(),
    ...) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object, layer = layer)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = "k")
    features <- list()
    for (i in as.numeric(
      x = names(x = table(object@kmeans.obj[[1]]$cluster))
    )) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning(
          "The following features are not present in the object: ",
          paste(missing.features, collapse = ", "),
          ifelse(
            test = search,
            yes = ", attempting to find updated synonyms",
            no = ", not searching for symbol synonyms"
          ),
          call. = FALSE,
          immediate. = TRUE
        )
        if (search) {
          tryCatch(
            expr = {
              updated.features <- UpdateSymbolList(
                symbols = missing.features,
                ...
              )
              names(x = updated.features) <- missing.features
              for (miss in names(x = updated.features)) {
                index <- which(x == miss)
                x[index] <- updated.features[miss]
              }
            },
            error = function(...) {
              warning(
                "Could not reach HGNC's gene names database",
                call. = FALSE,
                immediate. = TRUE
              )
            }
          )
          missing.features <- setdiff(x = x, y = rownames(x = object))
          if (length(x = missing.features) > 0) {
            warning(
              "The following features are still not present in the object: ",
              paste(missing.features, collapse = ", "),
              call. = FALSE,
              immediate. = TRUE
            )
          }
        }
      }
      return(intersect(x = x, y = rownames(x = object)))
    })
    cluster.length <- length(x = features)
  }
  if (!all(LengthCheck(values = features))) {
    warning(paste(
      "Could not find enough features in the object from the following feature lists:",
      paste(names(x = which(x = !LengthCheck(values = features)))),
      "Attempting to match case..."
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste(
      "The following feature lists do not have enough features present in the object:",
      paste(names(x = which(x = !LengthCheck(values = features)))),
      "exiting..."
    ))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- rowMeans(x = assay.data[pool, , drop = FALSE])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(
    x = data.avg + rnorm(n = length(data.avg)) / 1e+30,
    n = nbin,
    labels = FALSE,
    right = FALSE
  )
  names(x = data.cut) <- names(x = data.avg)

  scores <- bplapply(
    1:cluster.length,
    function(i) {
      features.use <- features[[i]]
      ctrl.use <- unlist(
        lapply(
          1:length(features.use),
          function(j) {
            data.cut[which(data.cut == data.cut[features.use[j]])]
          }
        )
      )
      ctrl.use <- names(
        sample(
          ctrl.use,
          size = min(ctrl * length(features.use), length(ctrl.use)),
          replace = FALSE
        )
      )
      ctrl.scores_i <- colMeans(x = assay.data[ctrl.use, , drop = FALSE])
      features.scores_i <- colMeans(
        x = assay.data[features.use, , drop = FALSE]
      )
      return(list(ctrl.scores_i, features.scores_i))
    },
    BPPARAM = BPPARAM
  )
  ctrl.scores <- do.call(rbind, lapply(scores, function(x) x[[1]]))
  features.scores <- do.call(rbind, lapply(scores, function(x) x[[2]]))

  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}
