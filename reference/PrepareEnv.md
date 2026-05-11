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
  ...
)
```

## Arguments

- envname:

  The name of the conda-compatible Python environment. If `NULL`, the
  environment name will be set to `"scop_env"`. Default is `NULL`.

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
  `"3.11-1"` on Windows.

- force:

  Whether to force recreation of the environment. If `TRUE`, the
  existing environment will be removed and recreated. Default is
  `FALSE`.

- modules:

  Optional Python dependency modules to install in addition to the
  default scientific stack. Supported values are `"scanpy"`, `"scvi"`,
  `"glue"`, `"scanorama"`, `"bbknn"`, `"celltypist"`, `"cellphonedb"`,
  `"magic"`, `"scrublet"`, `"sccoda"`, `"doubletdetection"`,
  `"doublet"`, `"palantir"`, `"scvelo"`, `"cellrank"`, `"wot"`,
  `"phate"`, `"pacmap"`, `"trimap"`, `"multimap"`, and `"scomm"`. If
  `NULL` or omitted in `PrepareEnv()`, the default environment is
  installed. The default excludes `"sccoda"` and `"scomm"` because their
  TensorFlow stacks are not compatible with the default JAX/scVI stack
  in the same environment; request them explicitly for scCODA/scOMM
  workflows.

- pip_options:

  Additional command line arguments to be passed to `uv`/`pip` when
  installing pip packages.

- ...:

  Additional arguments passed to package installation functions.
