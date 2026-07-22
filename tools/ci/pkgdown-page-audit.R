args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0L) "check" else args[[1L]]
full_build <- length(args) >= 2L && identical(args[[2L]], "true")
if (!mode %in% c("outputs", "check")) {
  stop("Usage: pkgdown-page-audit.R [outputs|check] [true|false]", call. = FALSE)
}

config <- yaml::read_yaml("_pkgdown.yml")
contents_or_empty <- function(section) {
  if (is.null(section$contents)) character() else section$contents
}

reference_topics <- unique(unlist(
  lapply(config$reference, contents_or_empty),
  use.names = FALSE
))
reference_topics <- reference_topics[
  is.character(reference_topics) & nzchar(reference_topics)
]

rd_files <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
rd_aliases <- lapply(rd_files, function(path) {
  lines <- readLines(path, warn = FALSE)
  entries <- grep("^\\\\(name|alias)\\{", lines, value = TRUE)
  sub("^\\\\(name|alias)\\{([^}]*)\\}.*$", "\\2", entries)
})
alias_to_page <- unlist(Map(
  function(aliases, path) {
    stats::setNames(
      rep(tools::file_path_sans_ext(basename(path)), length(aliases)),
      aliases
    )
  },
  rd_aliases,
  rd_files
))
missing_rd <- setdiff(reference_topics, names(alias_to_page))
if (length(missing_rd) > 0L) {
  stop(
    "Configured reference topics have no Rd name or alias: ",
    paste(missing_rd, collapse = ", "),
    call. = FALSE
  )
}
reference_pages <- unname(alias_to_page[reference_topics])

article_topics <- unique(unlist(
  lapply(config$articles, contents_or_empty),
  use.names = FALSE
))
article_topics <- article_topics[
  is.character(article_topics) & nzchar(article_topics)
]

expected_reference <- file.path("reference", paste0(reference_pages, ".html"))
expected_articles <- file.path("articles", paste0(article_topics, ".html"))

published <- tryCatch(
  system2(
    "git",
    c("ls-tree", "-r", "--name-only", "upstream/gh-pages"),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) character()
)
if (length(published) == 0L) {
  published <- tryCatch(
    system2(
      "git",
      c("ls-tree", "-r", "--name-only", "origin/gh-pages"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )
}
published <- gsub("\\\\", "/", published)

missing_reference <- reference_topics[!expected_reference %in% published]
missing_articles <- article_topics[!expected_articles %in% published]

if (identical(mode, "outputs")) {
  output_file <- Sys.getenv("GITHUB_OUTPUT", unset = "")
  output <- c(
    paste0("reference_topics=", paste(missing_reference, collapse = ",")),
    paste0("article_topics=", paste(missing_articles, collapse = ","))
  )
  if (nzchar(output_file)) {
    write(output, file = output_file, append = TRUE)
  } else {
    writeLines(output)
  }
  message(
    "Missing published pkgdown pages: ", length(missing_reference),
    " reference, ", length(missing_articles), " article"
  )
  quit(status = 0L)
}

local <- if (dir.exists("docs")) {
  list.files("docs", recursive = TRUE, all.files = TRUE, no.. = TRUE)
} else {
  character()
}
local <- gsub("\\\\", "/", local)
available <- if (full_build) local else unique(c(published, local))
missing <- c(
  expected_reference[!expected_reference %in% available],
  expected_articles[!expected_articles %in% available]
)
if (length(missing) > 0L) {
  stop(
    "Configured pkgdown pages are absent from gh-pages and this build:\n- ",
    paste(missing, collapse = "\n- "),
    call. = FALSE
  )
}
message(
  "Pkgdown page audit passed against ",
  if (full_build) "this full build" else "published and local output",
  ": ", length(expected_reference), " reference and ",
  length(expected_articles), " article pages"
)
