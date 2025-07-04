#' Try to evaluate an expression a set number of times before failing
#'
#' The function is used as a fail-safe if your R code sometimes works and sometimes
#' doesn't, usually because it depends on a resource that may be temporarily
#' unavailable. It tries to evaluate the expression `max_tries` times. If all the
#' attempts fail, it throws an error; if not, the evaluated expression is returned.
#'
#' @param expr The expression to be evaluated.
#' @param max_tries The maximum number of attempts to evaluate the expression before giving up. Default is set to 5.
#' @param error_message a string, additional custom error message you would like to be displayed when an error occurs.
#' @param retry_message a string, a message displayed when a new try to evaluate the expression would be attempted.
#'
#' @return This function returns the evaluated expression if successful,
#' otherwise it throws an error if all attempts are unsuccessful.
#' @export
#'
#' @examples
#' f <- function() {
#'   value <- runif(1, min = 0, max = 1)
#'   if (value > 0.5) {
#'     log_message("value is larger than 0.5")
#'     return(value)
#'   } else {
#'     log_message(
#'       "value is smaller than 0.5",
#'       message_type = "error"
#'     )
#'   }
#' }
#' f_evaluated <- try_get(expr = f())
#' print(f_evaluated)
try_get <- function(
    expr,
    max_tries = 5,
    error_message = "",
    retry_message = "Retrying...") {
  out <- simpleError("start")
  ntry <- 0
  while (inherits(out, "error")) {
    ntry <- ntry + 1
    out <- tryCatch(
      expr = eval.parent(substitute(expr)),
      error = function(error) {
        log_message(error)
        log_message("")
        log_message(error_message)
        Sys.sleep(1)
        return(error)
      }
    )
    if (inherits(out, "error") && ntry >= max_tries) {
      log_message(
        out,
        message_type = "error"
      )
    } else {
      if (!inherits(out, "error")) {
        break
      } else {
        log_message(retry_message)
      }
    }
  }
  return(out)
}

#' Download File from the Internet
#'
#' @md
#' @inheritParams utils::download.file
#' @param methods Methods to be used for downloading files.
#' The default is to try different download methods in turn until the download is successfully completed.
#' @param max_tries Number of tries for each download method.
#' @param verbose Logical value, default is `TRUE`.
#' Whether to print progress messages.
#' @param use_curl Logical value, default is `TRUE`.
#' Whether to use [curl::curl_download] to download the file.
#' @param ... Other arguments passed to [utils::download.file]
#'
#' @export
download <- function(
    url,
    destfile,
    methods = c("auto", "wget", "libcurl", "curl", "wininet", "internal"),
    max_tries = 2,
    use_curl = TRUE,
    verbose = TRUE,
    ...) {
  if (missing(url) || missing(destfile)) {
    log_message(
      "{.arg url} and {.arg destfile} must be both provided.",
      message_type = "error"
    )
  }
  ntry <- 0
  status <- NULL

  if (isTRUE(use_curl)) {
    tryCatch(
      {
        log_message(
          "Attempting download with {.val curl::curl_download}...",
          message_type = "info"
        )
        check_r("curl")

        h <- curl::new_handle()
        curl::handle_setopt(h, ssl_verifyhost = 0, ssl_verifypeer = 0)
        curl::handle_setheaders(h,
          "User-Agent" = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36",
          "Referer" = "https://guolab.wchscu.cn/AnimalTFDB4/",
          "Accept" = "*/*"
        )

        curl::curl_download(
          url = url,
          destfile = destfile,
          handle = h,
          quiet = !verbose
        )

        if (file.exists(destfile)) {
          content <- readLines(destfile, n = 1, warn = FALSE)
          if (!grepl("^<!doctype|^<html", tolower(content))) {
            log_message(
              "Download completed successfully using {.val curl::curl_download}",
              message_type = "success"
            )
            return(invisible(NULL))
          } else {
            unlink(destfile)
            log_message(
              "curl downloaded an HTML page, not data. Trying other methods.",
              "warning"
            )
          }
        }
      },
      error = function(e) {
        log_message(
          paste0("curl download failed: ", conditionMessage(e)),
          message_type = "warning"
        )
      }
    )
  }

  while (is.null(status)) {
    for (method in methods) {
      status <- tryCatch(
        expr = {
          suppressWarnings(
            utils::download.file(
              url = url,
              destfile = destfile,
              method = method,
              quiet = !verbose,
              ...
            )
          )
          status <- 1
        },
        error = function(error) {
          log_message(
            error,
            message_type = "warning"
          )
          log_message(
            "Cannot download from the url: ", url,
            message_type = "warning"
          )
          log_message(
            "Failed to download using {.val method}. Retry...",
            message_type = "warning"
          )
          Sys.sleep(1)
          return(NULL)
        }
      )
      if (!is.null(status)) {
        break
      }
    }
    ntry <- ntry + 1
    if (is.null(status) && ntry >= max_tries) {
      log_message(
        "Download failed.",
        message_type = "error"
      )
    }
  }
  return(invisible(NULL))
}

kegg_get <- function(url) {
  temp <- tempfile()
  on.exit(unlink(temp))
  download(url = url, destfile = temp)
  content <- as.data.frame(
    do.call(
      rbind,
      strsplit(readLines(temp), split = "\t")
    )
  )
  return(content)
}

rescale <- function(
    x,
    from = range(x, na.rm = TRUE, finite = TRUE),
    to = c(0, 1)) {
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  } else {
    return((x - from[1]) / diff(from) * diff(to) + to[1])
  }
}

zero_range <- function(
    x,
    tol = 1000 * .Machine$double.eps) {
  if (length(x) == 1) {
    return(TRUE)
  }
  if (length(x) != 2) {
    log_message(
      "x must be length 1 or 2",
      message_type = "error"
    )
  }
  if (any(is.na(x))) {
    return(NA)
  }
  if (x[1] == x[2]) {
    return(TRUE)
  }
  if (all(is.infinite(x))) {
    return(FALSE)
  }
  m <- min(abs(x))
  if (m == 0) {
    return(FALSE)
  }
  abs((x[1] - x[2]) / m) < tol
}

col2hex <- function(cname) {
  colMat <- grDevices::col2rgb(cname)
  grDevices::rgb(
    red = colMat[1, ] / 255,
    green = colMat[2, ] / 255,
    blue = colMat[3, ] / 255
  )
}

#' Invoke a function with a list of arguments
#' @param .fn A function, or function name as a string.
#' @param .args A list of arguments.
#' @param ... Other arguments passed to the function.
#' @param .env Environment in which to evaluate the call.
#' This will be most useful if .fn is a string, or the function has side-effects.
#'
#' @export
#'
#' @examples
#' f <- function(x, y) {
#'   x + y
#' }
#' invoke_fun(f, list(x = 1, y = 2))
#' invoke_fun("f", list(x = 1, y = 2))
#' invoke_fun("f", x = 1, y = 2)
invoke_fun <- function(
    .fn,
    .args = list(),
    ...,
    .env = rlang::caller_env()) {
  args <- c(.args, list(...))
  .bury <- c(".fn", "")
  if (rlang::is_null(.bury) || !length(args)) {
    if (rlang::is_scalar_character(.fn)) {
      .fn <- rlang::env_get(.env, .fn, inherit = TRUE)
    }
    call <- rlang::call2(.fn, !!!args)
    return(
      rlang::eval_bare(call, .env)
    )
  }
  if (!rlang::is_character(.bury, 2L)) {
    rlang::abort("`.bury` must be a character vector of length 2")
  }
  arg_prefix <- .bury[[2]]
  fn_nm <- .bury[[1]]
  buried_nms <- paste0(arg_prefix, seq_along(args))
  buried_args <- rlang::set_names(args, buried_nms)
  .env <- rlang::env(.env, !!!buried_args)
  args <- rlang::set_names(buried_nms, names(args))
  args <- rlang::syms(args)
  if (rlang::is_function(.fn)) {
    rlang::env_bind(.env, `:=`(!!fn_nm, .fn))
    .fn <- fn_nm
  }
  call <- rlang::call2(.fn, !!!args)
  return(
    rlang::eval_bare(call, .env)
  )
}

#' Implement similar functions to the \code{unnest} function in the tidyr package
#' @param data A data frame.
#' @param cols Columns to unnest.
#' @param keep_empty By default, you get one row of output for each element of the list your unchopping/unnesting.
#' This means that if there's a size-0 element (like \code{NULL} or an empty data frame),
#' that entire row will be dropped from the output.
#' If you want to preserve all rows,
#' use \code{keep_empty = TRUE} to replace size-0 elements with a single row of missing values.
#' @export
unnest_fun <- function(
    data,
    cols,
    keep_empty = FALSE) {
  if (nrow(data) == 0 || length(cols) == 0) {
    return(data)
  }
  for (col in cols) {
    col_expand <- unlist(data[[col]])
    expand_times <- sapply(data[[col]], length)
    if (isTRUE(keep_empty)) {
      data[[col]][expand_times == 0] <- NA
      col_expand <- unlist(data[[col]])
      expand_times[expand_times == 0] <- 1
    }
    data <- data[rep(seq_len(nrow(data)), times = expand_times), ]
    data[, col] <- col_expand
  }
  rownames(data) <- NULL
  return(data)
}

#' Capitalizes the characters
#' Making the first letter uppercase
#'
#' @param x A vector of character strings to be capitalized.
#' @param force_tolower Whether to force the remaining letters to be lowercase.
#' @export
#'
#' @examples
#' x <- c(
#'   "dna methylation",
#'   "rRNA processing",
#'   "post-Transcriptional gene silencing"
#' )
#' capitalize(x)
capitalize <- function(x, force_tolower = FALSE) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  if (!inherits(x, "character")) {
    log_message(
      "x must be the type of character.",
      message_type = "error"
    )
  }
  if (isTRUE(force_tolower)) {
    x <- paste(
      toupper(substr(x, 1, 1)),
      tolower(substr(x, 2, nchar(x))),
      sep = ""
    )
  } else {
    first_word <- sapply(strsplit(x, "\\s|-"), function(s) s[1])
    index <- which(first_word == tolower(first_word))
    x[index] <- paste(
      toupper(substr(x[index], 1, 1)),
      substr(x[index], 2, nchar(x[index])),
      sep = ""
    )
  }
  return(x)
}

str_wrap <- function(x, width = 80) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  x_wrap <- unlist(
    lapply(
      x,
      function(i) {
        paste0(strwrap(i, width = width), collapse = "\n")
      }
    )
  )
  return(x_wrap)
}

select_cells <- function(obj, celltypes, group.by) {
  metadata <- obj@meta.data
  cells_c <- c()
  for (celltype in celltypes) {
    cells_c <- c(
      cells_c,
      rownames(metadata[metadata[[group.by]] == celltype, ])
    )
  }
  return(cells_c)
}
