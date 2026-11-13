PIDC <- function(expMat, regulators = NULL, targets = NULL, logScale = FALSE, ncores = 1, diag = c("auto",
    "zero", "one"), verbose = interactive()) {
    if (verbose)
        thisutils::log_message("[1] Filtering and Normalizing data...")
    checkPIDCArgs(expMat = expMat, regulators = regulators, targets = targets, logScale = logScale, ncores = ncores)
    diag <- match.arg(diag)
    expMat <- as.matrix(expMat)
    if (logScale) {
        expMat <- log2(expMat + 1)
    }
    if (is.null(regulators)) {
        regulators <- row.names(expMat)
    }
    if (is.null(targets)) {
        targets <- row.names(expMat)
    }
    if (verbose)
        thisutils::log_message("[2] Discretizing expression matrix into bins...")
    expMat <- as.list(as.data.frame(t(expMat), check.names = FALSE))
    discret_list <- pbapply::pblapply(X = expMat, FUN = discretizeGene)
    regulators <- intersect(regulators, names(discret_list))
    targets <- intersect(targets, names(discret_list))
    if (verbose) {
        thisutils::log_message("[3] Calculating proportional unique contribution(PUC) ", "matrix using ",
            length(regulators), " regulators and ", length(targets), " targets...\n(This step may taken dozens ",
            "of hours, please be patient)")
    }
    if (ncores <= 1) {
        PUC_list <- pbapply::pblapply(X = targets, FUN = getPUC, regulators = regulators, discret_list = discret_list)
    }
    else {
        if (.Platform$OS.type == "unix") {
            cl <- parallel::makeCluster(spec = getOption("mc.cores", ncores), type = "FORK")
        }
        else {
            cl <- parallel::makeCluster(spec = getOption("mc.cores", ncores), type = "PSOCK")
        }
        parallel::clusterEvalQ(cl = cl, library(purrr))
        parallel::clusterEvalQ(cl = cl, library(dplyr))
        parallel::clusterEvalQ(cl = cl, library(pbapply))
        parallel::clusterExport(cl = cl, varlist = c("getPUC", "freqTable", "puc_per_target", "specific.information",
            "MI", "H"), envir = environment())
        PUC_list <- pbapply::pblapply(X = targets, FUN = getPUC, regulators = regulators, discret_list = discret_list,
            cl = cl)
        parallel::stopCluster(cl)
    }
    names(PUC_list) <- targets
    PUC_mat <- do.call(cbind, PUC_list)
    if (verbose)
        thisutils::log_message("[4] Summaring Uxy and calculating the weight matrix...")
    out_res <- FUxy(Uxy_mat = PUC_mat)
    if (diag == "auto") {
        return(out_res)
    }
    else if (diag == "zero") {
        diag(out_res) <- 0
    }
    else if (diag == "one") {
        diag(out_res) <- 1
    }
    out_res
}

checkPIDCArgs <- function(expMat, regulators, targets, logScale, ncores) {
    if (!is.matrix(expMat) && !is.array(expMat) && !is(expMat, "dgCMatrix")) {
        stop("The expMat must be a two-dimensional matrix where the row corresponds", " to a gene and each column corresponds to a sample")
    }
    if (length(dim(expMat)) != 2) {
        stop("The expMat must be a two-dimensional matrix where the row corresponds", " to a gene and each column corresponds to a sample")
    }
    if (is.null(rownames(expMat))) {
        stop("The expMat must contain the names of the genes as rownames.")
    }
    expMat <- as.matrix(expMat)
    if (!is.numeric(expMat))
        stop("The expMat contain non-numeric values.")
    if (any(is.na(expMat)))
        stop("The expMat contain NA.")
    gene_names <- unique(rownames(expMat))
    regulators <- unique(regulators)
    if (!is.null(regulators)) {
        if (!is.character(regulators))
            stop("The regulators must characters!")
        if (length(intersect(gene_names, regulators)) < 2) {
            stop("At least two regulators requiered, please check the names of the ", "expMat and the gene ids (names) in your regulators!")
        }
        else if (length(intersect(gene_names, regulators)) < length(regulators)) {
            ratio <- length(setdiff(regulators, gene_names))/length(regulators)
            warning(length(setdiff(regulators, gene_names)), " out of ", length(regulators), " (", round(ratio *
                100, 2), "%) regulators not found in your expMat!")
        }
    }
    targets <- unique(targets)
    if (!is.null(targets)) {
        if (!is.character(targets))
            stop("The targets must characters!")
        if (length(intersect(gene_names, targets)) < 1) {
            stop("At least one targets requiered, please check the names of the ", "expMat and the gene ids (names) in your targets!")
        }
        else if (length(intersect(gene_names, targets)) < length(targets)) {
            ratio <- length(setdiff(targets, gene_names))/length(targets)
            warning(length(setdiff(targets, gene_names)), " out of ", length(targets), " (", round(ratio *
                100, 2), "%) targets not found in your expMat!")
        }
    }
    if (!is.logical(logScale))
        stop("The logScale must be TRUE(T) or FALSE(F)!")
    if (!is.numeric(ncores) || ncores < 1)
        stop("The ncores must be positive numbers!")
}

discretizeGene <- function(x, method = c("uniform_width", "uniform_frequency")) {
    method <- match.arg(method)
    if (method == "uniform_width") {
        numBins <- floor(sqrt(length(x)))
        r <- range(x)
        b <- seq(from = r[1], to = r[2], length.out = (numBins + 1))
        X <- cut(x, breaks = b, include.lowest = TRUE)
    }
    else if ("uniform_frequency") {
        numBins <- floor(sqrt(length(x)))
        nrepl <- floor(length(x)/numBins)
        nplus <- sample(seq_len(numBins), length(x) - nrepl * numBins)
        nrep <- rep(nrepl, numBins)
        nrep[nplus] <- nrepl + 1
        X <- x[order(x)] <- rep(seq.int(numBins), nrep)
    }
    return(X)
}

freqTable <- function(x, y, method = c("ML", "Jeffreys", "Laplace", "SG", "minimax")) {
    method <- match.arg(method)
    if (is.list(x) || is.list(y)) {
        x <- unlist(x)
        y <- unlist(y)
    }
    tab <- table(x, y)
    if (method == "ML") {
        probs <- tab/sum(tab)
    }
    else if (method == "Jeffreys") {
        probs <- (tab + 0.5)/sum(tab + 0.5)
    }
    else if (method == "Laplace") {
        probs <- (tab + 1)/sum(tab + 1)
    }
    else if (method == "SG") {
        probs <- (tab + length(tab))/sum(tab + length(tab))
    }
    else if (method == "minimax") {
        probs <- (tab + sqrt(sum(tab))/length(tab))/sum(tab + sqrt(sum(tab))/length(tab))
    }
    probs
}

getPUC <- function(target, discret_list, regulators = names(discret_list)) {
    if (length(regulators) < 2)
        stop("At least two regulators needed!")
    if (length(target) != 1)
        stop("Only support one target each time!")
    if (is.null(names(discret_list)))
        stop("The names of discret_list must be gene symbols!")
    Z <- discret_list[target]
    Xi <- discret_list[regulators]
    p_XiZ <- purrr::map2(.x = Xi, .y = Z, .f = freqTable)
    p_Xi <- purrr::map(.x = p_XiZ, .f = rowSums)
    p_Z <- purrr::map(.x = p_XiZ, .f = colSums)
    I_XiZ <- purrr::pmap(.l = list(freqs = p_XiZ, p_i = p_Xi, p_z = p_Z), .f = MI)
    Ispec_XiZ <- purrr::pmap(.l = list(p_iz = p_XiZ, p_i = p_Xi, p_z = p_Z), .f = specific.information)
    U_xy <- puc_per_target(I_XiZ = I_XiZ, Ispec_XiZ = Ispec_XiZ, p_Z = p_Z)
    return(U_xy)
}

puc_per_target <- function(I_XiZ, Ispec_XiZ, p_Z) {
    gene_names <- names(Ispec_XiZ)
    Ispec_XiZ_bins <- lapply(purrr::transpose(Ispec_XiZ), unlist)
    Ispec_XiZ_bins_pos <- which(sapply(Ispec_XiZ_bins, sum) > 0)
    p_Z_bins <- lapply(purrr::transpose(p_Z), unlist)
    p_Z_bins_pos <- which(sapply(p_Z_bins, sum) > 0)
    if (!identical(Ispec_XiZ_bins_pos, p_Z_bins_pos)) {
        warning("The bins with values may not match between Ispec and p_Z!")
    }
    overlap_bins_pos <- union(Ispec_XiZ_bins_pos, p_Z_bins_pos)
    Ispec_XiZ_bins_filter <- Ispec_XiZ_bins[overlap_bins_pos]
    p_Z_bins_filter <- p_Z_bins[overlap_bins_pos]
    per_bin_redundancy <- purrr::map2(.x = Ispec_XiZ_bins_filter, .y = p_Z_bins_filter, .f = function(x,
        y) {
        if (is.null(names(x)) || is.null(names(y))) {
            stop("The Ispec and p_Z in each bin must have the same gene names!")
        }
        x_sort <- sort(x)
        n_gene <- length(x_sort)
        x_cumsum <- cumsum(x_sort)
        x_over_num <- n_gene - order(x_sort) - 1
        names(x_over_num) <- names(x_sort)
        x_over_sum <- x_sort * x_over_num
        x_value <- x_cumsum + x_over_sum
        y_value <- y[names(x_value)]
        out <- x_value * y_value
        out[names(x)]
    })
    check_gnames <- sapply(per_bin_redundancy, function(x, gname) {
        identical(names(x), gname)
    }, gname = gene_names)
    if (!all(check_gnames)) {
        per_bin_redundancy <- lapply(per_bin_redundancy, function(x, gnames) {
            x[gnames]
        }, gname = gene_names)
    }
    per_gene_redundancy <- sapply(purrr::transpose(per_bin_redundancy), function(x) {
        sum(unlist(x))
    })
    per_gene_puc <- length(gene_names) - 1 - per_gene_redundancy/unlist(I_XiZ)[names(per_gene_redundancy)]
    return(per_gene_puc)
}

specific.information <- function(p_iz, p_i, p_z, unit = c("log", "log2", "log10")) {
    unit <- match.arg(unit)
    p_i_z <- t(t(p_iz)/p_z)
    p_z_i <- t(p_iz/p_i)
    tmp <- t((log(1/p_z)) - (log(1/p_z_i))) * p_i_z
    if (unit == "log2")
        tmp <- t((log(1/p_z, 2)) - (log(1/p_z_i, 2))) * p_i_z
    if (unit == "log10")
        tmp <- t((log(1/p_z, 10)) - (log(1/p_z_i, 10))) * p_i_z
    colSums(tmp, na.rm = TRUE)
}

MI <- function(freqs, p_i, p_z, unit = c("log", "log2", "log10")) {
    if (length(dim(freqs)) != 2)
        stop("The dims of freqs must be 2!")
    unit <- match.arg(unit)
    MI <- H(p_i, unit = unit) + H(p_z, unit = unit) - H(freqs, unit = unit)
    MI
}

H <- function(freqs, unit = c("log", "log2", "log10")) {
    unit <- match.arg(unit)
    H <- -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
    if (unit == "log2")
        H <- H/log(2)
    if (unit == "log10")
        H <- H/log(10)
    H
}

FUxy <- function(Uxy_mat) {
    mat <- as.matrix(Uxy_mat)
    if (identical(row.names(mat), colnames(mat))) {
        mat <- mat + t(mat)
    }
    else if (nrow(mat) != ncol(mat)) {
        overlap_genes <- intersect(row.names(mat), colnames(mat))
        if (length(overlap_genes) > 0) {
            mat[overlap_genes, overlap_genes] <- (mat[overlap_genes, overlap_genes] + t(mat[overlap_genes,
                overlap_genes])) * 0.5
        }
    }
    else {
        stop("The row names and column names of mat is not identical!")
    }
    F_x <- t(apply(mat, 1, function(x) {
        (stats::ecdf(x))(x)
    }))
    colnames(F_x) <- colnames(mat)
    F_y <- apply(mat, 2, function(x) {
        (stats::ecdf(x))(x)
    })
    row.names(F_y) <- row.names(mat)
    F_xy <- (F_x + F_y) * 0.5
    return(F_xy)
}
