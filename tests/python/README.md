# PAGA compatibility regression

The default Euclidean neighbor search switches to PyNNDescent at 8192 cells.
`test_paga_neighbors.py` uses 8200 synthetic cells so it exercises this path,
then checks PAGA connectivity, layout, and graph embedding. It does not replace
Scanpy, PyNNDescent, or numerical results; only the R-specific logger is stubbed.
Install the dependencies listed in `.github/workflows/paga-python.yaml`, then run:

```sh
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
