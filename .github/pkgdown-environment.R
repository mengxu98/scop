plan <- jsonlite::read_json("pkgdown-plan.json", simplifyVector = TRUE)
installed <- installed.packages(fields = "RemoteSha")
installed <- installed[!duplicated(installed[, "Package"]), , drop = FALSE]
versions <- stats::setNames(
  paste(installed[, "Version"], ifelse(is.na(installed[, "RemoteSha"]), "", installed[, "RemoteSha"]), sep = ":"),
  installed[, "Package"]
)
previous <- if (!plan$full) {
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
