split_names <- function(value) {
  if (!nzchar(value)) character() else strsplit(value, ",", fixed = TRUE)[[1L]]
}

functions <- new.env(parent = emptyenv())
for (path in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  for (expression in parse(path)) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && is.call(expression[[3L]]) &&
      identical(expression[[3L]][[1L]], as.name("function"))) {
      functions[[as.character(expression[[2L]])]] <- expression[[3L]]
    }
  }
}

packages <- character()
visited <- character()
strings <- function(expression) {
  if (is.character(expression)) {
    return(expression)
  }
  if (is.call(expression) && identical(expression[[1L]], as.name("c"))) {
    return(unlist(lapply(as.list(expression)[-1L], strings), use.names = FALSE))
  }
  character()
}
visit <- function(expression, expand = TRUE) {
  if (is.expression(expression) || is.pairlist(expression)) {
    for (index in seq_along(expression)) {
      if (!identical(expression[[index]], quote(expr = ))) visit(expression[[index]], expand)
    }
    return(invisible(NULL))
  }
  if (!is.call(expression)) {
    return(invisible(NULL))
  }
  call <- expression[[1L]]
  if (identical(call, as.name("function"))) {
    visit(expression[[3L]], expand)
    return(invisible(NULL))
  }
  if (identical(call, as.name("if")) &&
    identical(expression[[2L]], quote(interactive()))) {
    if (length(expression) == 4L) visit(expression[[4L]], expand)
    return(invisible(NULL))
  }
  if (is.symbol(call)) {
    name <- as.character(call)
    if (name %in% c("::", ":::")) packages <<- union(packages, as.character(expression[[2L]]))
    if (name %in% c("check_r", "get_namespace_fun", "getExportedValue", "packageDescription") && length(expression) > 1L) {
      packages <<- union(packages, strings(expression[[2L]]))
    }
    if (name %in% c("library", "require", "requireNamespace") && length(expression) > 1L) {
      package <- expression[[2L]]
      if (is.symbol(package) || is.character(package)) packages <<- union(packages, as.character(package))
    }
    if (expand && exists(name, functions, inherits = FALSE) && !name %in% visited) {
      visited <<- c(visited, name)
      methods <- ls(functions)
      methods <- methods[methods == name | startsWith(methods, paste0(name, "."))]
      for (method in methods) visit(functions[[method]], expand = FALSE)
    }
  }
  for (index in seq_along(expression)) {
    if (!identical(expression[[index]], quote(expr = )) && !is.symbol(expression[[index]])) {
      visit(expression[[index]], expand)
    }
  }
  invisible(NULL)
}

for (topic in split_names(Sys.getenv("PKGDOWN_TOPICS"))) {
  example <- tempfile(fileext = ".R")
  tools::Rd2ex(file.path("man", paste0(topic, ".Rd")), out = example)
  if (file.exists(example)) visit(parse(example))
  unlink(example)
}
for (article in split_names(Sys.getenv("PKGDOWN_ARTICLES"))) {
  lines <- readLines(file.path("vignettes", paste0(article, ".Rmd")), warn = FALSE)
  active <- FALSE
  code <- character()
  for (line in lines) {
    if (grepl("^```\\{r([ ,}]|$)", line)) {
      active <- !grepl("eval\\s*=\\s*(FALSE|F)([, }]|$)", line)
    } else if (grepl("^```\\s*$", line)) {
      active <- FALSE
    } else if (active) {
      code <- c(code, line)
    }
  }
  visit(parse(text = code))
}

sources <- c(
  thisplot = "github::mengxu98/thisplot",
  thisutils = "github::mengxu98/thisutils",
  harmony = "github::immunogenomics/harmony",
  monocle3 = "github::cole-trapnell-lab/monocle3",
  SeuratData = "github::satijalab/seurat-data",
  spacexr = "github::dmcable/spacexr",
  CytoTRACE2 = "github::digitalcytometry/cytotrace2/cytotrace2_r"
)
if ("jinworks/SpatialCellChat" %in% packages || "SpatialCellChat" %in% packages) {
  packages <- union(packages, c("KlugerLab/ALRA", "JEFworks-Lab/MERINGUE"))
}
github <- packages[grepl("/", packages, fixed = TRUE) & !grepl("::", packages, fixed = TRUE)]
sources[sub(".*/", "", github)] <- paste0("github::", github)
base_packages <- rownames(installed.packages(priority = "base"))
packages <- sort(setdiff(packages, c(base_packages, "scop")))
references <- vapply(packages, function(package) {
  if (package %in% names(sources)) {
    return(sources[[package]])
  }
  if (grepl("/", package, fixed = TRUE) && !grepl("::", package, fixed = TRUE)) {
    return(paste0("github::", package))
  }
  package
}, "", USE.NAMES = FALSE)
execute <- identical(Sys.getenv("PKGDOWN_EXECUTE"), "true")
required <- if (execute) "deps::.,any::sessioninfo" else "any::sessioninfo"
extra <- unique(c("any::pkgdown", if (execute) c("local::.", "any::testthat"), references))
output <- c(paste0("packages=", required), paste0("extra=", paste(extra, collapse = ",")))
writeLines(output)
destination <- Sys.getenv("GITHUB_OUTPUT")
if (nzchar(destination)) write(output, file = destination, append = TRUE)
