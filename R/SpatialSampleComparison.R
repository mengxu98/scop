#' Summarize spatial composition or neighborhoods for each sample
#'
#' Composition counts labelled observations, not deconvolved cells. Neighborhood
#' summaries average the observed fractions across target cells within each
#' sample, giving every target cell equal weight. Empty neighborhoods stay NA.
#' Missing categories in composition are genuine zero counts. No cell-level
#' hypothesis test is performed.
#'
#' @param object A Seurat object or cell metadata data.frame with sample design
#'   columns and unique cell IDs as row names. A metadata table needs no artificial
#'   expression matrix, including for imaging cytometry data.
#' @param group.by Cell type or domain metadata column for composition.
#' @param sample.by Sample or section identifier. Required.
#' @param condition.by Experimental condition. Each sample has exactly one.
#' @param subject.by Independent animal/patient identifier. NULL asserts that
#'   every sample is an independent subject; supply this for repeated sections.
#' @param profile Optional output of SpatialNeighborhoodProfile. Use its cells
#'   argument to select targets while retaining their full tissue context.
#' @return A data.frame with sample, condition, subject, group, lower, radius,
#'   estimate and n_observations. source and measure attributes document scope.
#' @seealso SpatialSampleComparison, SpatialSamplePlot
#' @export
SpatialSampleSummary <- function(object, group.by, sample.by, condition.by,
                                  subject.by = NULL, profile = NULL) {
  if (!inherits(object, "Seurat") && !is.data.frame(object)) stop("object must be a Seurat object or cell metadata table", call. = FALSE)
  for (nm in list(group.by, sample.by, condition.by)) validate_scalar_string(nm, "design column")
  if (!is.null(subject.by)) validate_scalar_string(subject.by, "subject.by")
  meta <- if (is.data.frame(object)) object else object[[]]
  if (!nrow(meta) || is.null(rownames(meta)) || anyDuplicated(rownames(meta))) stop("Cell metadata must have unique row IDs", call. = FALSE)
  columns <- c(group.by, sample.by, condition.by, subject.by)
  if (!all(columns %in% names(meta))) stop("Design columns are missing", call. = FALSE)
  if (any(vapply(meta[, columns, drop = FALSE], function(x) anyNA(x) || any(!nzchar(as.character(x))), logical(1)))) {
    stop("Design and group labels must be nonmissing and nonempty", call. = FALSE)
  }
  design <- data.frame(sample = as.character(meta[[sample.by]]),
    condition = as.character(meta[[condition.by]]),
    subject = as.character(meta[[subject.by %||% sample.by]]), stringsAsFactors = FALSE)
  design <- unique(design)
  if (anyDuplicated(design$sample)) stop("Each sample must map to exactly one condition and subject", call. = FALSE)
  if (is.null(profile)) {
    tab <- as.data.frame(table(sample = as.character(meta[[sample.by]]),
      group = as.character(meta[[group.by]])), stringsAsFactors = FALSE)
    names(tab)[3] <- "count"
    tab$n_observations <- stats::ave(tab$count, tab$sample, FUN = sum)
    tab$estimate <- tab$count / tab$n_observations
    tab$lower <- tab$radius <- NA_real_
    measure <- "fraction of labelled observations"
    source <- list(cells = rownames(meta), group.by = group.by)
  } else {
    required <- c("cell_id", "sample", "group", "lower", "radius", "fraction", "total")
    if (!is.data.frame(profile) || !nrow(profile) || !all(required %in% names(profile))) {
      stop("profile must be a nonempty SpatialNeighborhoodProfile table", call. = FALSE)
    }
    idx <- match(profile$cell_id, rownames(meta))
    if (anyNA(idx) || anyNA(profile$sample) ||
        !identical(as.character(profile$sample), as.character(meta[[sample.by]][idx]))) {
      stop("Profile cell IDs/sample labels do not match the sample design", call. = FALSE)
    }
    if (!is.numeric(profile$fraction) || any(!is.na(profile$fraction) &
        (!is.finite(profile$fraction) | profile$fraction < 0 | profile$fraction > 1)) ||
        any(!is.finite(profile$lower) | !is.finite(profile$radius)) ||
        any(profile$radius <= profile$lower) || anyNA(profile$group)) {
      stop("Profile fractions or radius intervals are invalid", call. = FALSE)
    }
    keys <- c("cell_id", "group", "lower", "radius")
    if (anyDuplicated(profile[, keys])) stop("Profile contains duplicate target/group/radius rows", call. = FALSE)
    strata <- unique(profile[, c("group", "lower", "radius"), drop = FALSE])
    tab <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i) {
      ss <- strata[i, ]
      rows <- profile[profile$group == ss$group & profile$lower == ss$lower & profile$radius == ss$radius, ]
      do.call(rbind, lapply(design$sample, function(sample) {
        values <- rows$fraction[rows$sample == sample]
        data.frame(sample = sample, group = as.character(ss$group), lower = ss$lower,
          radius = ss$radius, count = NA_real_,
          estimate = if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE),
          n_observations = sum(!is.na(values)), stringsAsFactors = FALSE)
      }))
    }))
    measure <- "mean neighborhood fraction per target observation"
    source <- list(profile_source = attr(profile, "source"),
      profile_parameters = attr(profile, "parameters"), target_cells = unique(profile$cell_id))
  }
  out <- merge(design, tab, by = "sample", sort = FALSE)
  out <- out[, c("sample", "subject", "condition", "group", "lower", "radius", "count", "n_observations", "estimate")]
  attr(out, "source") <- c(source, list(sample.by = sample.by,
    condition.by = condition.by, subject.by = subject.by,
    independent_samples_asserted = is.null(subject.by)))
  attr(out, "measure") <- measure
  out
}

#' Compare spatial summaries using independent subjects
#'
#' Average repeated sections equally within each subject and condition before
#' inference. Independent contrasts use Welch's t test; paired contrasts use a
#' paired t test on complete subjects. Effects are treatment minus reference in
#' the original fraction units. Insufficient replication or undefined standard
#' error produces not_tested with a reason, never a cell-level replacement.
#' BH correction is applied across all testable rows in this call.
#'
#' @param summary Output of SpatialSampleSummary.
#' @param contrast Explicit c(reference, treatment) condition labels.
#' @param paired Whether subjects are measured in both conditions. Subjects in
#'   both conditions are rejected for an independent contrast.
#' @param conf.level Confidence level for the mean difference.
#' @return A list with comparisons, subject_values, sample_values, parameters
#'   and source. comparisons includes effects, intervals, p/FDR, replication
#'   counts, excluded unpaired subjects and test status. This is an unadjusted
#'   two-condition analysis, not a covariate-adjusted or causal model.
#' @export
SpatialSampleComparison <- function(summary, contrast, paired = FALSE, conf.level = 0.95) {
  required <- c("sample", "subject", "condition", "group", "lower", "radius", "estimate")
  if (!is.data.frame(summary) || !nrow(summary) || !all(required %in% names(summary))) {
    stop("summary must be a nonempty SpatialSampleSummary table", call. = FALSE)
  }
  validate_scalar_flag(paired, "paired")
  if (!is.character(contrast) || length(contrast) != 2L || anyNA(contrast) ||
      any(!nzchar(contrast)) || anyDuplicated(contrast)) stop("contrast must be c(reference, treatment)", call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1L || !is.finite(conf.level) ||
      conf.level <= 0 || conf.level >= 1) stop("conf.level must be between zero and one", call. = FALSE)
  if (!all(contrast %in% summary$condition)) stop("Both contrast conditions must be present", call. = FALSE)
  if (!is.numeric(summary$estimate) || any(!is.na(summary$estimate) & !is.finite(summary$estimate))) {
    stop("estimate must be numeric and finite or NA", call. = FALSE)
  }
  for (nm in c("sample", "subject", "condition", "group")) {
    if (anyNA(summary[[nm]]) || any(!nzchar(as.character(summary[[nm]])))) stop("Summary design labels are invalid", call. = FALSE)
  }
  if (anyDuplicated(summary[, c("sample", "group", "lower", "radius")])) {
    stop("Summary has duplicate sample/group/radius rows", call. = FALSE)
  }
  mapping <- unique(summary[, c("sample", "subject", "condition")])
  if (anyDuplicated(mapping$sample)) stop("Samples have inconsistent subject/condition mappings", call. = FALSE)
  dat <- summary[summary$condition %in% contrast, , drop = FALSE]
  if (!paired && length(intersect(dat$subject[dat$condition == contrast[1]], dat$subject[dat$condition == contrast[2]]))) {
    stop("Subjects occur in both conditions; use paired = TRUE", call. = FALSE)
  }
  strata <- unique(dat[, c("group", "lower", "radius"), drop = FALSE])
  subject_values <- comparisons <- vector("list", nrow(strata))
  equal <- function(x, y) (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
  for (i in seq_len(nrow(strata))) {
    ss <- strata[i, ]
    rows <- dat[dat$group == ss$group & equal(dat$lower, ss$lower) & equal(dat$radius, ss$radius), ]
    units <- unique(rows[, c("subject", "condition"), drop = FALSE])
    units$estimate <- vapply(seq_len(nrow(units)), function(j) {
      v <- rows$estimate[rows$subject == units$subject[j] & rows$condition == units$condition[j]]
      if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    }, numeric(1))
    units$n_samples <- vapply(seq_len(nrow(units)), function(j) {
      sum(rows$subject == units$subject[j] & rows$condition == units$condition[j] & !is.na(rows$estimate))
    }, integer(1))
    x <- units[units$condition == contrast[2] & !is.na(units$estimate), ]
    y <- units[units$condition == contrast[1] & !is.na(units$estimate), ]
    excluded <- 0L
    if (paired) {
      common <- intersect(x$subject, y$subject)
      excluded <- length(union(x$subject, y$subject)) - length(common)
      x <- x[match(common, x$subject), , drop = FALSE]
      y <- y[match(common, y$subject), , drop = FALSE]
    }
    units$included_in_contrast <- !is.na(units$estimate) & units$subject %in% union(x$subject, y$subject)
    subject_values[[i]] <- cbind(units, ss[rep(1L, nrow(units)), , drop = FALSE])
    effect <- if (nrow(x) && nrow(y)) mean(x$estimate) - mean(y$estimate) else NA_real_
    test <- NULL
    reason <- if (min(nrow(x), nrow(y)) < 2L) "fewer than two independent subjects per condition or complete pairs" else NA_character_
    if (is.na(reason)) {
      test <- tryCatch(stats::t.test(x$estimate, y$estimate, paired = paired, conf.level = conf.level),
        error = function(e) e)
      if (inherits(test, "error")) { reason <- conditionMessage(test); test <- NULL }
    }
    comparisons[[i]] <- cbind(ss, data.frame(reference = contrast[1], treatment = contrast[2],
      effect = effect, conf_low = if (is.null(test)) NA_real_ else test$conf.int[1],
      conf_high = if (is.null(test)) NA_real_ else test$conf.int[2],
      p_value = if (is.null(test)) NA_real_ else test$p.value,
      n_reference = nrow(y), n_treatment = nrow(x), n_unpaired_excluded = excluded,
      status = if (is.null(test)) "not_tested" else "tested", reason = reason, stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, comparisons)
  out$q_value <- stats::p.adjust(out$p_value, method = "BH")
  list(comparisons = out, subject_values = do.call(rbind, subject_values), sample_values = summary,
    parameters = list(method = if (paired) "paired_t" else "welch_t", contrast = contrast,
      paired = paired, conf.level = conf.level, aggregation = "equal section means within subject and condition",
      correction = "BH across all testable rows", measure = attr(summary, "measure")),
    source = attr(summary, "source"))
}

#' Plot the independent subject values used in a spatial comparison
#'
#' @param result Output of SpatialSampleComparison.
#' @return A ggplot showing subject values; lines connect complete pairs when
#'   paired = TRUE. No backend is rerun and no significance is inferred in plotting.
#' @export
SpatialSamplePlot <- function(result) {
  if (!is.list(result) || !is.data.frame(result$subject_values)) stop("result must be a SpatialSampleComparison result", call. = FALSE)
  dat <- result$subject_values
  if ("included_in_contrast" %in% names(dat)) dat <- dat[dat$included_in_contrast, , drop = FALSE]
  dat$condition <- factor(dat$condition, levels = result$parameters$contrast)
  dat$panel <- ifelse(is.na(dat$radius), as.character(dat$group),
    paste(dat$group, paste0("(", dat$lower, ", ", dat$radius, "]")))
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$condition, y = .data$estimate))
  if (isTRUE(result$parameters$paired)) p <- p + ggplot2::geom_line(ggplot2::aes(group = .data$subject), alpha = 0.35, na.rm = TRUE)
  p + ggplot2::geom_point(ggplot2::aes(color = .data$condition), size = 2, na.rm = TRUE) +
    ggplot2::facet_wrap(~panel, scales = "free_y") + ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = result$parameters$measure %||% "Fraction", color = "Condition")
}
