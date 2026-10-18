arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments)) match.arg(arguments, c("dependencies", "environment", "build")) else "build"

if (mode == "dependencies") {
  split_names <- function(value) {
    if (!nzchar(value)) character() else strsplit(value, ",", fixed = TRUE)[[1L]]
  }

  functions <- new.env(parent = emptyenv())
  exports <- grep("^export\\(", readLines("NAMESPACE", warn = FALSE), value = TRUE)
  exports <- gsub('"', "", sub("^export\\((.*)\\)$", "\\1", exports))
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
      if ((expand || !name %in% exports) && exists(name, functions, inherits = FALSE) && !name %in% visited) {
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
    # Temporary runtime pin until CRAN includes r-lib/pak#917. Without the fix,
    # pak can unload processx's SIGCHLD handler before parallel restores it.
    pak = "github::r-lib/pak@6f4c90a69a6dafcc5b1fbfd715feba22055e77df",
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
  extra <- unique(c("any::pkgdown", if (execute) c("local::.", sources[["pak"]]), references))
  output <- c(paste0("packages=", required), paste0("extra=", paste(extra, collapse = ",")))
  writeLines(output)
  destination <- Sys.getenv("GITHUB_OUTPUT")
  if (nzchar(destination)) write(output, file = destination, append = TRUE)
  quit(status = 0L)
}

installed_versions <- function() {
  installed <- installed.packages(fields = "RemoteSha")
  versions <- stats::setNames(
    paste(installed[, "Version"], ifelse(is.na(installed[, "RemoteSha"]), "", installed[, "RemoteSha"]), sep = ":"),
    installed[, "Package"]
  )
  versions[!duplicated(names(versions))]
}
plan <- jsonlite::read_json("pkgdown-plan.json", simplifyVector = TRUE)
if (mode == "environment") {
  versions <- installed_versions()
  previous <- if (!plan$full) {
    if (system2("git", c("fetch", "--no-tags", "--depth=1", "origin", plan$published)) != 0L) {
      stop("Unable to fetch the published website revision")
    }
    jsonlite::fromJSON(system2("git", c("show", paste0(plan$published, ":.pkgdown-build.json")), stdout = TRUE))
  } else {
    NULL
  }
  changed <- FALSE
  if (!is.null(previous)) {
    common <- intersect(names(versions), names(previous$packages))
    changed <- !identical(previous$R, as.character(getRversion())) ||
      is.null(previous$packages) || any(versions[common] != unlist(previous$packages[common]))
  }
  if (changed) {
    plan$full <- TRUE
    plan$execute <- TRUE
    plan$topics <- tools::file_path_sans_ext(list.files("man", pattern = "[.]Rd$"))
    articles <- list.files("vignettes", pattern = "[.]Rmd$", recursive = TRUE)
    plan$articles <- tools::file_path_sans_ext(articles[!grepl("(^|/)_", articles)])
    plan$home <- plan$news <- plan$reference_index <- plan$articles_index <- TRUE
    plan$reasons <- c(plan$reasons, "R or installed dependency versions changed")
    for (field in c("topics", "articles", "deleted", "reasons")) plan[[field]] <- I(plan[[field]])
    jsonlite::write_json(plan, "pkgdown-plan.json", auto_unbox = TRUE, pretty = TRUE)
    write(c(
      paste0("PKGDOWN_TOPICS=", paste(plan$topics, collapse = ",")),
      paste0("PKGDOWN_ARTICLES=", paste(plan$articles, collapse = ",")),
      "PKGDOWN_EXECUTE=true"
    ), file = Sys.getenv("GITHUB_ENV"), append = TRUE)
  }
  write(paste0("changed=", tolower(changed)), file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)
  quit(status = 0L)
}

check_pages <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing generated pages: ", paste(missing, collapse = ", "))
  errors <- lapply(paths, function(path) {
    page <- xml2::read_html(path)
    messages <- xml2::xml_text(xml2::xml_find_all(page, "//*[contains(concat(' ', normalize-space(@class), ' '), ' r-err ')]"))
    if (length(messages)) paste(path, paste(unique(messages), collapse = "\n"), sep = ":\n") else character()
  })
  errors <- unlist(errors, use.names = FALSE)
  if (length(errors)) stop("Example execution failed:\n", paste(errors, collapse = "\n\n"))
}

timings <- data.frame(page = character(), seconds = numeric())
record <- function(name, expression) {
  start <- proc.time()[["elapsed"]]
  on.exit({
    timings <<- rbind(timings, data.frame(page = name, seconds = proc.time()[["elapsed"]] - start))
    utils::write.csv(timings, "pkgdown-timings.csv", row.names = FALSE)
  })
  force(expression)
}
capture.output(sessioninfo::session_info(), file = "pkgdown-session.txt")
installed <- installed.packages(fields = "RemoteSha")
installed <- installed[!duplicated(installed[, "Package"]), , drop = FALSE]
utils::write.table(installed[, c("Package", "Version", "RemoteSha")],
  file = "pkgdown-session.txt", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE
)
options(Ncpus = 2L, future.globals.maxSize = 4 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential", gc = TRUE)
pkg <- pkgdown::as_pkgdown(".")
build <- function() {
  if (plan$full) record("initialise", pkgdown::init_site(pkg))
  if (plan$home) record("home", pkgdown::build_home(pkg))
  if (plan$news) record("news", pkgdown::build_news(pkg, preview = FALSE))
  if (plan$reference_index) record("reference/index", pkgdown::build_reference_index(pkg))
  for (topic in plan$topics) {
    record(paste0("reference/", topic), {
      pkgdown::build_reference(pkg, topics = topic, examples = TRUE, run_dont_run = FALSE, devel = FALSE)
      check_pages(file.path(pkg$dst_path, "reference", paste0(topic, ".html")))
    })
  }
  if (plan$articles_index) record("articles/index", pkgdown::build_articles_index(pkg))
  for (article in plan$articles) {
    record(paste0("articles/", article), {
      pkgdown::build_article(article, pkg = pkg)
      check_pages(file.path(pkg$dst_path, "articles", paste0(article, ".html")))
    })
  }
  if (plan$full) {
    record("site-metadata", pkgdown::build_site(
      pkg,
      lazy = TRUE, examples = TRUE, run_dont_run = FALSE, install = FALSE,
      new_process = FALSE, preview = FALSE
    ))
  } else {
    record("search", pkgdown::build_search(pkg))
  }
}
record(if (plan$full) "full-site" else "incremental-site", build())
if (length(plan$deleted) && any(file.exists(file.path("docs", plan$deleted)))) {
  stop("Deleted documentation pages remain in website output")
}
file.create("docs/.nojekyll")
versions <- installed_versions()
previous <- if (file.exists("docs/.pkgdown-build.json")) {
  jsonlite::read_json("docs/.pkgdown-build.json")$packages
} else {
  list()
}
previous[names(versions)] <- as.list(versions)
jsonlite::write_json(
  list(source = plan$source, R = as.character(getRversion()), packages = previous),
  "docs/.pkgdown-build.json",
  auto_unbox = TRUE, pretty = TRUE
)
summary <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summary)) {
  write(
    c(
      "", "| Page | Seconds |", "| --- | ---: |",
      sprintf("| %s | %.2f |", timings$page, timings$seconds)
    ),
    file = summary, append = TRUE
  )
}
