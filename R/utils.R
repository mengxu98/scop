#' Download File from the Internet
#'
#' @md
#' @inheritParams utils::download.file
#' @param methods Methods to be used for downloading files.
#' The default is to try different download methods in turn until the download is successfully completed.
#' @param use_httr Logical value, default is `FALSE`.
#' Whether to use the `httr` package to download files.
#' @param verbose Logical value, default is `TRUE`.
#' @param ... Other arguments passed to [utils::download.file]
#' @param max_tries Number of tries for each download method.
#'
#' @export
download <- function(
    url,
    destfile,
    methods = c("auto", "wget", "libcurl", "curl", "wininet", "internal"),
    use_httr = FALSE,
    verbose = TRUE,
    ...,
    max_tries = 2) {
  if (missing(url) || missing(destfile)) {
    log_message(
      "{.arg url} and {.arg destfile} must be both provided.",
      message_type = "error"
    )
  }

  ntry <- 0
  status <- NULL

  if (isTRUE(use_httr)) {
    tryCatch(
      {
        check_r("httr")

        log_message(
          "Attempting download with {.val httr::GET}...",
          message_type = "info"
        )

        resp <- httr::GET(
          url,
          httr::add_headers(
            "User-Agent" = "Mozilla/5.0",
            "Referer" = dirname(url),
            "Accept" = "*/*"
          ),
          httr::write_disk(destfile, overwrite = TRUE)
        )

        if (httr::status_code(resp) == 200) {
          log_message(
            "Download completed successfully using {.val httr::GET}",
            message_type = "success"
          )
          return(invisible(NULL))
        } else {
          log_message(
            paste0("HTTP error: ", httr::status_code(resp), " - ", httr::http_status(resp)$message),
            message_type = "warning"
          )
        }
      },
      error = function(e) {
        log_message(
          paste0("httr::GET failed: ", conditionMessage(e)),
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
            paste0("Cannot download from the url: ", url),
            message_type = "warning"
          )
          log_message(
            paste0("Failed to download using ", method, ". Retry..."),
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

col2hex <- function(cname) {
  colMat <- grDevices::col2rgb(cname)
  grDevices::rgb(
    red = colMat[1, ] / 255,
    green = colMat[2, ] / 255,
    blue = colMat[3, ] / 255
  )
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
