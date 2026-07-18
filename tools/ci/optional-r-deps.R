collect_check_r_packages <- function(root = ".") {
  collect_strings <- function(expr) {
    if (is.character(expr)) return(expr)
    if (is.call(expr) && identical(as.character(expr[[1]]), "c")) {
      return(unlist(lapply(as.list(expr)[-1], collect_strings), use.names = FALSE))
    }
    character()
  }

  walk <- function(expr) {
    if (!is.call(expr)) return(character())
    packages <- if (identical(as.character(expr[[1]]), "check_r") && length(expr) >= 2L) {
      collect_strings(expr[[2]])
    } else {
      character()
    }
    unique(c(packages, unlist(lapply(as.list(expr)[-1], walk), use.names = FALSE)))
  }

  source_files <- list.files(
    file.path(root, "R"),
    pattern = "[.]R$",
    full.names = TRUE
  )
  expressions <- unlist(lapply(source_files, function(path) {
    as.list(parse(path, keep.source = FALSE))
  }), recursive = FALSE)
  unique(unlist(lapply(expressions, walk), use.names = FALSE))
}

optional_r_packages <- function(root = ".") {
  description <- read.dcf(file.path(root, "DESCRIPTION"))
  fields <- intersect(
    c("Depends", "Imports", "Suggests", "LinkingTo", "Remotes"),
    colnames(description)
  )
  values <- description[1, fields]
  values <- values[!is.na(values)]
  description_packages <- trimws(unlist(strsplit(values, ",", fixed = TRUE)))
  description_packages <- sub("\\s*\\(.*$", "", description_packages)
  description_packages <- sub("^(github|gitlab)::", "", description_packages)
  description_packages <- description_packages[
    nzchar(description_packages) & description_packages != "R"
  ]

  unique(c(description_packages, collect_check_r_packages(root)))
}

manifest_path <- Sys.getenv("SCOP_DEP_MANIFEST")
if (nzchar(manifest_path)) {
  writeLines(sort(optional_r_packages()), manifest_path)
}
