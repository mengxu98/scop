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
#' @param max_tries Number of tries for each download method.
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
