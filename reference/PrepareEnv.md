# Prepare the python environment

Prepare the python environment by installing the required dependencies
and setting up the environment.

## Usage

``` r
PrepareEnv(
  envname = NULL,
  conda = "auto",
  miniconda_repo = "https://repo.anaconda.com/miniconda",
  version = if (is_windows()) "3.11-1" else "3.10-1",
  force = FALSE,
  modules = NULL,
  pip_options = character(),
  components = "python",
  cores = NULL,
  verbose = TRUE,
  configure_runtime = TRUE,
  ...
)
```

## Arguments

- envname:

  Conda-compatible Python environment. If `NULL`, use the current
  environment selection. `PrepareEnv()` selects a separate module
  profile when required.

- conda:

  The path or command name of a conda-compatible executable (`conda`,
  `mamba`, or `micromamba`). Use `"auto"` to allow automatically finding
  an appropriate environment manager. If `"micromamba"` is requested and
  micromamba is not available on `PATH`, a package-managed micromamba is
  downloaded automatically.

- miniconda_repo:

  Repository URL for miniconda. Default is
  <https://repo.anaconda.com/miniconda>.

- version:

  The Python version. Default is `"3.10-1"` on macOS and Unix and
  `"3.11-1"` on Windows. Existing compatible environments are reused;
  trajectory package pins are selected for the actual Python
  3.10/3.11/3.12 profile.

- force:

  Whether to force recreation of the environment. If `TRUE`, the
  existing environment will be removed and recreated.

- modules:

  Optional Python dependency modules to install in addition to the
  default scientific stack. If `NULL` or omitted in `PrepareEnv()`, the
  default environment is installed. The default excludes `"phate"`,
  `"magic"`, `"cell2location"`, `"cell2fate"`, `"commot"`, `"sccoda"`,
  and `"scomm"` because these workflows use long-running or incompatible
  dependency stacks. Cell2fate and COMMOT are prepared in standalone
  `"cell2fate_env"` and `"commot_env"` environments with Python 3.9.
  SCENIC+ uses standalone `"scenicplus_env"` or the smaller
  `"scenicplus_io_env"` for `.h5mu` imports. `"scenic"` is also excluded
  from the default environment and is prepared in `"scenic_env"` by
  default because SCENIC requires an older Python/numpy stack. PHATE and
  MAGIC are prepared explicitly in `"scprep_env"` with pandas 2.0.3,
  separately from the default scientific stack. scCODA and scOMM use
  `"tensorflow_env"` to keep TensorFlow dependencies separate from JAX.
  `"scmalignantfinder"` is prepared in a standalone
  `"scmalignantfinder_env"` because its upstream package pins an older
  scanpy/anndata stack. On Windows, the default also excludes `"scvi"`,
  `"glue"`, and `"multimap"` because those upstream stacks are more
  reliable when requested explicitly for method-specific workflows.

- pip_options:

  Additional command line arguments to be passed to `uv`/`pip` when
  installing pip packages.

- components:

  Components to prepare. Supported values are `"python"`, `"r"`, and
  `"all"`; `"all"` expands to both Python and R. The default,
  `"python"`, preserves the existing Python-only behavior. The R
  component collects packages declared through `check_r()` and
  checks/installs them, including optional workflow dependencies not
  listed in `DESCRIPTION`. CHOIR is R-only, so
  `PrepareEnv(modules = "choir")` prepares its optional R backend
  without creating a Python environment.

- configure_runtime:

  Whether to select the prepared interpreter for embedded Python and
  configure the R process. Set to `FALSE` for subprocess-only use.

- cores:

  Number of workers for R-package installation. Use `NULL` (the default)
  to let pak select its worker count automatically.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to package installation functions.
