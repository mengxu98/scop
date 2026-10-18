# Exercise the real exit-time regression from https://github.com/r-lib/pak/pull/917.
# This must run in a fresh Rscript with the same default libraries as the website.
expected_sha <- "6f4c90a69a6dafcc5b1fbfd715feba22055e77df"
description <- utils::packageDescription("pak")
if (!identical(description[["RemoteSha"]], expected_sha)) {
  stop("The website runtime must use pak with the pinned native-exit fix.")
}
cat("Runtime pak:", as.character(utils::packageVersion("pak")), "\n")
cat("Runtime pak path:", find.package("pak"), "\n")
cat("Runtime pak revision:", description[["RemoteSha"]], "\n")

# Load parallel before creating pak's private processx finalizer. Its first fork
# saves processx's SIGCHLD handler, then its exit finalizer restores that handler.
loadNamespace("parallel")
pak:::load_private_package("processx", "c_")
process <- pak:::pkg_data$ns$processx$process$new("true")
process$wait(timeout = 10000)
stopifnot(!process$is_alive(), identical(process$get_exit_status(), 0L))
job <- parallel::mcparallel(1L)
result <- parallel::mccollect(job)
stopifnot(identical(unname(result), list(1L)))

# Do not unload DLLs or use quit(): a normal Rscript exit is the test assertion.
cat("Completed process/fork checks; testing normal Rscript shutdown.\n")
