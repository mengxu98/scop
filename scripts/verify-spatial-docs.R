pkgload::load_all(".", quiet = TRUE, compile = FALSE)
topics <- c("ReadSpatialData", "RunStandardWorkflow", "RunSpatialVariableFeatures")
stopifnot(all(topics %in% getNamespaceExports("scop")))
for (topic in topics) {
  path <- file.path("man", paste0(topic, ".Rd"))
  stopifnot(file.exists(path))
  tools::checkRd(tools::parse_Rd(path))
  examples <- tempfile(fileext = ".R")
  tools::Rd2ex(path, out = examples)
  if (file.exists(examples)) parse(examples)
  unlink(examples)
}
pkgdown::build_reference(topics = topics, examples = FALSE)
for (article in "spatial-platform-workflows") {
  pkgdown::build_article(article)
}
