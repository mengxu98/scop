# scop

# scop 0.1.9

* **Enhancements**: Update `CellScoring()` and `AddModuleScore2()` functions. Now, new parameters `cores` and `verbose` have been added. The `AddModuleScore2()` function no longer uses the `BiocParallel::bpparam()` function to enable parallelization, but [parallelize_fun](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html), and the `cores` parameter is used to control the number of cores in [parallelize_fun](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html).

# scop 0.1.5

* **Bug Fixes**: Fix error for `RunPalantir()` function, see [#23](https://github.com/mengxu98/scop/issues/23).

# scop 0.1.4

* **Enhancements**: Update `.onAttach()`, now `.onAttach()` will print more information about conda and Python.
* **Enhancements**: Update `PrepareEnv()` function for easy add or update a conda environments and install Python packages.
* **New Functionality**: Added `ListEnv()` and `RemoveEnv()` functions for easy management of conda environment and Python packages.

# scop 0.1.3

* **New Functionality**: Added `TACSPlot()` function for creating FACS-like plots. Please refer to [manuscript](https://doi.org/10.1016/j.immuni.2018.04.015) and [code](https://github.com/maehrlab/thymusatlastools2/blob/f8b51ad684d56b2eeda780787eb9ad4ff3003eef/R/data_handling_seurat.R#L271) for specific information.

# scop 0.0.9

* **Bug Fixes**: Corrected an issue in `PrepareEnv()` function. The default Python version is now set to `3.10`, which ensures that Python-dependent functions like `RunPAGA()` and `RunSCVELO()` function run correctly.

# scop 0.0.6

* **New Functionality**: Add `GetAssayData5()` function, a reimplementation of `GetAssayData()`, for compatibility with Seurat v5 `Assay` objects.
* **Enhancements**: Updated `GetFeaturesData()` and `AddFeaturesData()` function to support retrieving and adding feature metadata.

# scop 0.0.5

* **Data**: Updated the `pancreas_sub` and `panc8_sub` test datasets to the Seurat v5 object format.
