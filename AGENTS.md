# AGENTS.md

## Cursor Cloud specific instructions

`scop` is an R package for single-cell and spatial omics analysis. It is a
standard R package: R code in `R/`, a compiled C++/Rcpp backend in `src/`
(40+ translation units using RcppArmadillo/RcppEigen/RcppParallel + OpenMP),
bundled example datasets in `data/`, and `testthat` tests in `tests/testthat/`.
The authoritative CI pipeline is `.github/workflows/R-CMD-check.yaml`
(`R CMD check` on Linux/macOS/Windows with only *hard* dependencies).

### Environment baseline (already provisioned in the VM snapshot)

- R 4.6.1 is installed system-wide (`/usr/bin/R`), plus the compiler toolchain
  and system libraries needed by the C++ backend and dependencies (HDF5, GDAL,
  GEOS, GSL, GLPK, FFTW, fontconfig/harfbuzz/freetype, libxml2, etc.).
- `/usr/lib/R/etc/Rprofile.site` pins package repos to the Posit Public Package
  Manager **binary** mirrors (CRAN for Ubuntu noble + Bioconductor 3.23) and
  sets `Ncpus`. This is why installs are fast — do not remove it, or packages
  will fall back to slow source compilation. `pak` adds any remaining
  Bioconductor source repos itself.
- `/usr/local/lib/R/site-library` is owned by `ubuntu`, so package installs do
  **not** need `sudo`.
- All *hard* dependencies (DESCRIPTION `Depends`/`Imports`/`LinkingTo`, incl.
  Seurat, Signac, and Bioconductor packages), `pak`, `testthat`, `pkgload`, and
  the compiled `scop` package are already installed. `thisutils`/`thisplot`
  come from GitHub (`Remotes:`).
- The update script only refreshes hard dependencies (`pak::local_install_deps`).
  It does **not** reinstall `scop` — after editing R or `src/` code, rebuild
  manually (see below).

### Build / install the package

- `MAKEFLAGS="-j4" R CMD INSTALL --no-multiarch --with-keep.source .`
- `src/Makevars` is generated from `src/Makevars.in` by the `./configure` script
  at build time (detects RcppParallel flags); do not edit `src/Makevars` directly.
- A full rebuild recompiles the whole C++ backend (a few minutes). `pkgload`'s
  `load_all()` also recompiles, so prefer `R CMD INSTALL` when you only need the
  installed package.

### Run (the "application" is the R package)

- Quick end-to-end smoke test using the bundled dataset:
  ```r
  library(scop)
  data(pancreas_sub)
  pancreas_sub <- RunStandardWorkflow(pancreas_sub)   # normalize -> HVF -> PCA -> cluster -> UMAP
  CellDimPlot(pancreas_sub, group.by = c("CellType", "SubCellType"), reduction = "UMAP")
  ```
- Other bundled datasets live in `data/` (e.g. `panc8_sub`, `pbmcmultiome_sub`,
  `visium_*`). See `README.md` for the full workflow gallery.

### Test

- Canonical developer flow (mirrors how tests expect the package namespace to be
  loaded, and makes `local_mocked_bindings` work):
  ```r
  pkgload::load_all(".")
  testthat::test_local(".", filter = "cpp-backend-safety")   # omit filter to run all
  ```
- Many tests reference **non-exported** functions directly. When running a test
  file without `load_all`, pass `env = asNamespace("scop")` to
  `testthat::test_dir(...)`, otherwise internal C++ helpers (e.g.
  `plage_dense`, `as_matrix`) are "not found".
- Only *hard* dependencies are installed, so tests gated on optional `Suggests`
  packages (GSVA, CIBERSORT, BiocNeighbors, harmony, scran, ...) **skip**
  cleanly — that is expected, not a failure. Install a specific Suggests with
  `pak::pkg_install("<pkg>")` if you need to exercise those paths.

### Lint / check

- The project's "lint" is `R CMD check` (the `R-CMD-check` workflow). A full
  check needs the large `Suggests` set; for a fast static/lint pass without
  installing every optional package:
  ```sh
  R CMD build --no-build-vignettes --no-manual .
  _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --no-vignettes \
    --no-tests --no-examples scop_0.9.0.tar.gz
  ```
  This validates syntax, NAMESPACE, Rd docs, code/doc consistency, and compiled
  code (Status: OK).
