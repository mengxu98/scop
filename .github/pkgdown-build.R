arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments)) match.arg(arguments, c("dependencies", "python", "environment", "build")) else "build"

installed_versions <- function() {
  installed <- installed.packages(fields = "RemoteSha")
  versions <- stats::setNames(
    paste(installed[, "Version"], ifelse(is.na(installed[, "RemoteSha"]), "", installed[, "RemoteSha"]), sep = ":"),
    installed[, "Package"]
  )
  versions[!duplicated(names(versions))]
}
plan <- jsonlite::read_json("pkgdown-plan.json", simplifyVector = TRUE)

if (mode == "dependencies") {
  code <- expression()
  dormant <- character()
  topics <- unique(c(plan$topics, tools::file_path_sans_ext(list.files("man", pattern = "\\.Rd$"))))
  for (topic in topics) {
    path <- tempfile(fileext = ".R")
    rd <- tools::parse_Rd(file.path("man", paste0(topic, ".Rd")))
    tools::Rd2ex(rd, out = path)
    expressions <- if (file.exists(path)) parse(path) else expression()
    if (topic %in% plan$topics) code <- c(code, expressions)
    if (!length(expressions) && any(vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\examples"), logical(1)))) {
      aliases <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\alias"), rd)
      dormant <- c(dormant, vapply(aliases, function(x) paste(unlist(x), collapse = ""), ""))
    }
    unlink(path)
  }
  for (article in plan$articles) {
    lines <- readLines(file.path("vignettes", paste0(article, ".Rmd")), warn = FALSE)
    active <- FALSE
    chunk <- character()
    for (line in lines) {
      if (grepl("^```\\{r([ ,}]|$)", line)) {
        active <- !grepl("eval\\s*=\\s*(FALSE|F)([, }]|$)", line)
      } else if (grepl("^```\\s*$", line)) {
        active <- FALSE
      } else if (active) {
        chunk <- c(chunk, line)
      }
    }
    code <- c(code, parse(text = chunk))
  }
  python <- list()
  noninteractive_code <- function(x) {
    if (is.call(x) && identical(x[[1L]], as.name("if")) && identical(x[[2L]], quote(interactive()))) {
      return(if (length(x) == 4L) noninteractive_code(x[[4L]]) else quote(NULL))
    }
    if (is.call(x) && (identical(x[[1L]], as.name("PrepareEnv")) ||
        identical(x[[1L]], quote(scop::PrepareEnv)))) {
      components <- as.list(x)$components
      if (is.null(components) || "python" %in% eval(components, baseenv())) {
        python[[length(python) + 1L]] <<- x
      }
    }
    if (!is.call(x) && !is.expression(x) && !is.pairlist(x)) return(x)
    parts <- as.list(x)
    for (i in seq_along(parts)) {
      if (!identical(parts[[i]], quote(expr = ))) parts[i] <- list(noninteractive_code(parts[[i]]))
    }
    if (is.call(x)) as.call(parts) else if (is.pairlist(x)) as.pairlist(parts) else as.expression(parts)
  }
  code <- noninteractive_code(code)
  literals <- function(x) {
    if (missing(x)) return(character())
    if (is.character(x)) return(x)
    if (!is.call(x) && !is.expression(x) && !is.pairlist(x)) return(character())
    unlist(lapply(as.list(x), literals), use.names = FALSE)
  }
  selectors <- unique(literals(code))
  callbacks <- function(x) {
    if (missing(x) || (!is.call(x) && !is.expression(x) && !is.pairlist(x))) return(character())
    parts <- as.list(x)
    selected <- if (is.call(x) && identical(x[[1L]], quote(list))) {
      Filter(is.symbol, parts[names(parts) %in% selectors])
    }
    nested <- if (is.call(x)) parts[-1L] else parts
    unique(c(vapply(selected, as.character, ""), unlist(lapply(nested, callbacks), use.names = FALSE)))
  }
  python <- unique(python)
  manifest <- jsonlite::read_json("pkgdown-plan.json")
  manifest$python <- lapply(python, function(x) paste(deparse(x), collapse = "\n"))
  jsonlite::write_json(manifest, "pkgdown-plan.json", auto_unbox = TRUE, pretty = TRUE)
  path <- tempfile()
  saveRDS(python, path, version = 2)
  write(c(paste0("python=", tolower(length(python) > 0L)), paste0("python_key=", tools::md5sum(path))),
    file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)
  unlink(path)
  definitions <- list()
  for (path in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
    for (x in parse(path)) {
      if (is.call(x) && identical(x[[1L]], quote(`<-`)) && is.symbol(x[[2L]]) &&
          is.call(x[[3L]]) && identical(x[[3L]][[1L]], quote(`function`))) {
        definitions[[as.character(x[[2L]])]] <- x[[3L]][[3L]]
      }
    }
  }
  namespace <- parseNamespaceFile(".", ".")
  methods <- namespace$S3methods
  method_names <- ifelse(is.na(methods[, 3L]), paste(methods[, 1L], methods[, 2L], sep = "."), methods[, 3L])
  names(method_names) <- methods[, 1L]
  called <- intersect(all.names(code), names(definitions))
  dormant <- setdiff(intersect(dormant, namespace$exports), called)
  visited <- character()
  while (length(called)) {
    called <- unique(c(called, unname(method_names[names(method_names) %in% called])))
    called <- setdiff(intersect(called, names(definitions)), visited)
    if (!length(called)) break
    bodies <- as.expression(definitions[called])
    code <- c(code, bodies)
    visited <- c(visited, called)
    selected <- callbacks(bodies)
    dormant <- setdiff(dormant, selected)
    calls <- unlist(lapply(bodies, function(body) {
      fun <- function(...) NULL
      body(fun) <- body
      codetools::findGlobals(fun, merge = FALSE)$functions
    }), use.names = FALSE)
    called <- setdiff(intersect(c(calls, selected), names(definitions)), c(visited, dormant))
  }
  strings <- function(x) {
    if (is.character(x)) return(x)
    if (is.call(x) && identical(x[[1L]], quote(c))) return(unlist(lapply(as.list(x)[-1L], strings), use.names = FALSE))
    character()
  }
  declared <- function(x) {
    if (missing(x) || !is.call(x)) return(character())
    package <- if (is.symbol(x[[1L]]) && as.character(x[[1L]]) %in% c("check_r", "get_namespace_fun") && length(x) >= 2L) strings(x[[2L]]) else character()
    unique(c(package, unlist(lapply(as.list(x)[-1L], declared), use.names = FALSE)))
  }
  optional <- unique(unlist(lapply(code, declared), use.names = FALSE))
  path <- tempfile(fileext = ".R")
  writeLines(unlist(lapply(code, deparse), use.names = FALSE), path)
  dependencies <- renv::dependencies(path, root = dirname(path), progress = FALSE, errors = "fatal")
  unlink(path)
  packages <- unique(c(dependencies$Package, optional))
  packages <- packages[order(!grepl("/", packages, fixed = TRUE))]
  package_names <- sub("@.*$", "", sub("^.*/", "", packages))
  packages <- packages[!duplicated(package_names)]
  description <- read.dcf("DESCRIPTION")
  packages <- setdiff(packages, c(description[1L, "Package"], rownames(installed.packages(priority = "base"))))
  locked <- if (file.exists("renv.lock")) jsonlite::read_json("renv.lock")$Packages else list()
  packages <- vapply(sort(packages), function(package) {
    record <- locked[[package]]
    if (identical(record$Source, "GitHub")) {
      return(paste0("github::", paste(c(record$RemoteUsername, record$RemoteRepo, record$RemoteSubdir), collapse = "/")))
    }
    package
  }, "", USE.NAMES = FALSE)
  write(paste0("extra=", paste(unique(packages), collapse = ",")),
    file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)
  quit(status = 0L)
}

if (mode == "python") {
  prepared <- lapply(plan$python, function(code) {
    request <- parse(text = code)[[1L]]
    cache <- callr::r(function(request) {
      args <- as.list(match.call(scop::PrepareEnv, request))[-1L]
      args <- lapply(args, eval, envir = baseenv())
      if (is.null(args$conda) || identical(args$conda, "auto")) args$conda <- "micromamba"
      do.call(scop::PrepareEnv, args)
      getOption("scop_env_cache")
    }, list(request), show = TRUE, user_profile = FALSE)
    write(paste0("RETICULATE_CONDA=", cache$conda), file = Sys.getenv("GITHUB_ENV"), append = TRUE)
    cache
  })
  for (cache in prepared) {
    valid <- getFromNamespace("exist_python_pkgs", "scop")(
      cache$requirements$packages, envname = cache$envname, conda = cache$conda, verbose = FALSE
    )
    if (!length(valid) || !all(valid)) stop("Python environment requirements are not satisfied: ", paste(names(valid)[!valid], collapse = ", "))
  }
  quit(status = 0L)
}

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
  write(paste0("changed=", tolower(changed)), file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)
  quit(status = 0L)
}

check_pages <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing generated pages: ", paste(missing, collapse = ", "))
  errors <- lapply(paths, function(path) {
    page <- xml2::read_html(path)
    messages <- xml2::xml_text(xml2::xml_find_all(page, "//*[contains(concat(' ', normalize-space(@class), ' '), ' r-err ')]"))
    output <- xml2::xml_text(xml2::xml_find_all(page,
      "//*[contains(concat(' ', normalize-space(@class), ' '), ' r-msg ') or contains(concat(' ', normalize-space(@class), ' '), ' r-wrn ')]"
    ))
    messages <- c(messages, output[grepl("\\berror\\b.*\\bskip", output, ignore.case = TRUE)])
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
pkg <- pkgdown::as_pkgdown(".")
build <- function() {
  on.exit(parallel:::cleanup(shutdown = TRUE), add = TRUE)
  if (plan$full) {
    pkgdown::build_site(
      pkg, examples = TRUE, run_dont_run = FALSE, install = FALSE,
      new_process = FALSE, preview = FALSE, devel = FALSE
    )
  } else {
    if (plan$home) record("home", pkgdown::build_home(pkg))
    if (plan$news) record("news", pkgdown::build_news(pkg, preview = FALSE))
    if (length(plan$topics)) {
      record("reference", pkgdown::build_reference(
        pkg, topics = plan$topics, examples = TRUE, run_dont_run = FALSE, devel = FALSE
      ))
    } else if (plan$reference_index) {
      record("reference/index", pkgdown::build_reference_index(pkg))
    }
    if (plan$articles_index) record("articles/index", pkgdown::build_articles_index(pkg))
    for (article in plan$articles) {
      record(paste0("articles/", article), {
        pkgdown::build_article(article, pkg = pkg)
      })
    }
    record("search", pkgdown::build_search(pkg))
  }
  check_pages(c(
    file.path(pkg$dst_path, "reference", paste0(plan$topics, ".html", recycle0 = TRUE)),
    file.path(pkg$dst_path, "articles", paste0(plan$articles, ".html", recycle0 = TRUE))
  ))
}
record(if (plan$full) "full-site" else "incremental-site", build())
capture.output(sessioninfo::session_info(), file = "pkgdown-session.txt", append = TRUE)
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
      "", "| Stage | Seconds |", "| --- | ---: |",
      sprintf("| %s | %.2f |", timings$page, timings$seconds)
    ),
    file = summary, append = TRUE
  )
}
