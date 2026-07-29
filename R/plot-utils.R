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

resolve_legacy_plot_theme <- function(
  theme_use = "theme_scop",
  theme_args = list()
) {
  if (identical(theme_use, "theme_scop")) {
    theme_use <- "theme_this"
  }
  theme_fun <- tryCatch(
    get(theme_use, mode = "function"),
    error = function(e) NULL
  )
  if (is.null(theme_fun)) {
    return(ggplot2::theme_bw())
  }
  do.call(theme_fun, theme_args)
}

resolve_method_plot_theme <- function(
  theme_use = "theme_scop",
  theme_args = list()
) {
  if (is.null(theme_use)) {
    return(ggplot2::theme_minimal())
  }
  if (inherits(theme_use, "theme")) {
    return(theme_use)
  }
  theme_fun <- if (is.character(theme_use)) {
    get(theme_use, mode = "function", inherits = TRUE)
  } else {
    theme_use
  }
  do.call(theme_fun, theme_args)
}
