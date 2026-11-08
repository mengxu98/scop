# PAGA compatibility regression

The default Euclidean neighbor search switches to PyNNDescent at 8192 cells.
`test_paga_neighbors.py` uses 8200 synthetic cells so it exercises this path,
then checks PAGA connectivity, layout, and graph embedding. It does not replace
Scanpy, PyNNDescent, or numerical results; only the R-specific logger is stubbed.
The existing Ubuntu job in `.github/workflows/R-CMD-check.yaml` runs this
regression in an isolated Python 3.10 environment. No separate workflow or job
is required. To run locally, install the same dependencies and execute:

```sh
python -m pip install 'scanpy==1.11.3' 'scipy==1.15.3' 'pynndescent==0.6.0' 'igraph==0.11.9' 'numpy==1.26.4' 'pandas==2.2.0' 'anndata==0.11.4' 'zarr<3' 'numba==0.66.0' 'llvmlite==0.48.0' 'scikit-learn==1.7.0' 'matplotlib==3.10.8'
python -m pip check
NUMBA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests/python -v
```

For the Seurat roundtrip, install scop's R dependencies and use a fresh R session:

```sh
RETICULATE_PYTHON=/absolute/path/to/python SCOP_RUN_PAGA_PYTHON_TESTS=true \
NUMBA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
Rscript -e 'testthat::test_local(filter="paga-python-integration")'
```

The R integration test bypasses only `PrepareEnv()` provisioning, using the
preinstalled Python environment. Seurat conversion, Scanpy neighbors/PAGA,
layout and conversion back to Seurat all run normally. Environment-profile
pins are separately checked by `test-prepare-env-requirements.R`.
