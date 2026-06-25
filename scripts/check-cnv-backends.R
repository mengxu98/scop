#!/usr/bin/env Rscript

options(scop_env_init = FALSE)

flag <- tolower(Sys.getenv("SCOP_RUN_CNV_BACKENDS", "false"))
if (!flag %in% c("1", "true", "yes")) {
  message("Set SCOP_RUN_CNV_BACKENDS=true to run optional CNV backend smoke checks.")
  quit(status = 0)
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not installed: ", pkg, call. = FALSE)
  }
}

read_required_rds <- function(envvar, label) {
  path <- Sys.getenv(envvar, "")
  if (!nzchar(path) || !file.exists(path)) {
    stop("Set ", envvar, " to an RDS file for ", label, call. = FALSE)
  }
  readRDS(path)
}

require_package("scop")
require_package("Seurat")

srt <- read_required_rds("SCOP_CNV_SRT_RDS", "the input Seurat object")
gene_order <- if (nzchar(Sys.getenv("SCOP_CNV_GENE_ORDER_RDS", ""))) {
  readRDS(Sys.getenv("SCOP_CNV_GENE_ORDER_RDS"))
} else {
  NULL
}
reference_by <- Sys.getenv("SCOP_CNV_REFERENCE_BY", "")
reference <- Sys.getenv("SCOP_CNV_REFERENCE", "")
if (!nzchar(reference_by) || !nzchar(reference)) {
  stop("Set SCOP_CNV_REFERENCE_BY and SCOP_CNV_REFERENCE for reference-cell methods.", call. = FALSE)
}

run_one <- function(method, expr) {
  message("Running CNV backend: ", method)
  force(expr)
  message("Finished CNV backend: ", method)
}

if (requireNamespace("copykat", quietly = TRUE)) {
  run_one("copykat", scop::RunCNV(srt, method = "copykat", genome = "hg38", verbose = FALSE))
}
if (requireNamespace("fastCNV", quietly = TRUE)) {
  run_one("fastCNV", scop::RunCNV(
    srt,
    method = "fastCNV",
    reference.by = reference_by,
    reference = reference,
    genome = "hg38",
    verbose = FALSE
  ))
}
if (requireNamespace("SCEVAN", quietly = TRUE)) {
  run_one("scevan", scop::RunCNV(srt, method = "scevan", genome = "hg38", verbose = FALSE))
}
if (requireNamespace("infercnv", quietly = TRUE)) {
  run_one("infercnv", scop::RunCNV(
    srt,
    method = "infercnv",
    reference.by = reference_by,
    reference = reference,
    gene_order = gene_order,
    genome = "hg38",
    verbose = FALSE
  ))
}
if (requireNamespace("numbat", quietly = TRUE)) {
  allele_counts <- read_required_rds("SCOP_CNV_ALLELE_COUNTS_RDS", "Numbat allele counts")
  reference_counts <- read_required_rds("SCOP_CNV_REFERENCE_COUNTS_RDS", "Numbat reference counts")
  run_one("numbat", scop::RunCNV(
    srt,
    method = "numbat",
    allele_counts = allele_counts,
    reference_counts = reference_counts,
    genome = "hg38",
    verbose = FALSE
  ))
}
if (requireNamespace("CaSpER", quietly = TRUE)) {
  loh <- read_required_rds("SCOP_CNV_LOH_RDS", "CaSpER LOH signal")
  cytoband <- read_required_rds("SCOP_CNV_CYTOBAND_RDS", "CaSpER cytoband table")
  run_one("casper", scop::RunCNV(
    srt,
    method = "casper",
    reference.by = reference_by,
    reference = reference,
    gene_order = gene_order,
    loh = loh,
    cytoband = cytoband,
    genome = "hg38",
    verbose = FALSE
  ))
}

message("Optional CNV backend smoke checks completed.")
