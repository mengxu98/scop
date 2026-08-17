combine_plot_list <- function(
  plots,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
) {
  if (!isTRUE(combine)) {
    return(plots)
  }
  if (length(plots) > 1L) {
    return(patchwork::wrap_plots(
      plotlist = plots,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow
    ))
  }
  plots[[1L]]
}

matrix_to_long <- function(mat, row_name, col_name, value_name) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c(row_name, col_name, value_name)
  df[[row_name]] <- factor(df[[row_name]], levels = rownames(mat))
  df[[col_name]] <- factor(df[[col_name]], levels = colnames(mat))
  df
}

# "theme_scop" is a scop alias of thisplot::theme_this. thisplot::StatPlot()
# looks up character theme_use in the thisplot namespace, so the string
# "theme_scop" cannot be found there. Return the function instead.
resolve_plot_theme_use <- function(theme_use) {
  if (identical(theme_use, "theme_scop")) {
    return(thisplot::theme_this)
  }
  if (identical(theme_use, "theme_spatial")) {
    return(theme_spatial)
  }
  theme_use
}

# Build a ggplot2 theme from theme_use (name, function, or theme object).
# allow_null = TRUE returns NULL when theme_use is NULL instead of the fallback.
apply_plot_theme <- function(
  theme_use = "theme_scop",
  theme_args = list(),
  fallback = ggplot2::theme_bw,
  allow_null = FALSE
) {
  if (is.null(theme_use)) {
    if (isTRUE(allow_null)) {
      return(NULL)
    }
    return(do.call(fallback, list()))
  }
  if (inherits(theme_use, "theme")) {
    return(theme_use)
  }
  theme_use <- resolve_plot_theme_use(theme_use)
  theme_fun <- if (is.function(theme_use)) {
    theme_use
  } else if (is.character(theme_use) && length(theme_use) == 1L) {
    tryCatch(
      get(theme_use, mode = "function", inherits = TRUE),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  if (is.null(theme_fun)) {
    return(do.call(fallback, list()))
  }
  fmls <- names(formals(theme_fun))
  if (length(theme_args) > 0L && !is.null(names(theme_args))) {
    keep <- setdiff(fmls, "...")
    if ("..." %in% fmls) {
      keep <- union(keep, names(formals(ggplot2::theme)))
    }
    theme_args <- theme_args[names(theme_args) %in% keep]
  }
  do.call(theme_fun, theme_args)
}
