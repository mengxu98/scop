# Package index

### Single-Cell Omics analysis Pipeline

- [`scop`](https://mengxu98.github.io/scop/reference/scop-package.md)
  [`scop-package`](https://mengxu98.github.io/scop/reference/scop-package.md)
  : Single-Cell Omics analysis Pipeline
- [`scop_logo()`](https://mengxu98.github.io/scop/reference/scop_logo.md)
  : scop logo

### Package Management

- [`check_python()`](https://mengxu98.github.io/scop/reference/check_python.md)
  : Check and install python packages
- [`check_r()`](https://mengxu98.github.io/scop/reference/check_r.md) :
  Check and install R packages
- [`conda_install()`](https://mengxu98.github.io/scop/reference/conda_install.md)
  : Enhanced conda installation
- [`conda_python()`](https://mengxu98.github.io/scop/reference/conda_python.md)
  : Find the path to Python associated with a conda environment
- [`env_exist()`](https://mengxu98.github.io/scop/reference/env_exist.md)
  : Check if a conda environment exists
- [`env_info()`](https://mengxu98.github.io/scop/reference/env_info.md)
  : Print environment information
- [`env_requirements()`](https://mengxu98.github.io/scop/reference/env_requirements.md)
  : Python environment requirements
- [`exist_python_pkgs()`](https://mengxu98.github.io/scop/reference/exist_python_pkgs.md)
  : Check if the python package exists in the environment
- [`find_conda()`](https://mengxu98.github.io/scop/reference/find_conda.md)
  : Find an appropriate conda binary
- [`install_miniconda2()`](https://mengxu98.github.io/scop/reference/install_miniconda2.md)
  : Enhanced miniconda installation
- [`installed_python_pkgs()`](https://mengxu98.github.io/scop/reference/installed_python_pkgs.md)
  : Show all the python packages in the environment
- [`ListEnv()`](https://mengxu98.github.io/scop/reference/ListEnv.md) :
  List conda environments
- [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
  : Prepare the python environment
- [`remove_python()`](https://mengxu98.github.io/scop/reference/remove_python.md)
  : Remove Python packages from conda environment
- [`remove_r()`](https://mengxu98.github.io/scop/reference/remove_r.md)
  : Check and remove R packages
- [`RemoveEnv()`](https://mengxu98.github.io/scop/reference/RemoveEnv.md)
  : Remove a conda environment

### scop pipeline

- [`integration_scop()`](https://mengxu98.github.io/scop/reference/integration_scop.md)
  : The integration_scop function
- [`standard_scop()`](https://mengxu98.github.io/scop/reference/standard_scop.md)
  : Standard workflow for scop

### Data Processing

- [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md)
  : Convert an anndata object to a seurat object using reticulate

- [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md)
  : Cell scoring

- [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  :

  Check and report the type of data in `Seurat` object

- [`CheckDataList()`](https://mengxu98.github.io/scop/reference/CheckDataList.md)
  :

  Check and preprocess a list of `Seurat` objects

- [`CheckDataMerge()`](https://mengxu98.github.io/scop/reference/CheckDataMerge.md)
  : Check and preprocess a merged seurat object

- [`CycGenePrefetch()`](https://mengxu98.github.io/scop/reference/CycGenePrefetch.md)
  : Prefetch cell cycle genes

- [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
  : Find the default reduction name in a Seurat object

- [`FetchDataZero()`](https://mengxu98.github.io/scop/reference/FetchDataZero.md)
  : FetchData but with zeroes for unavailable genes

- [`GeneConvert()`](https://mengxu98.github.io/scop/reference/GeneConvert.md)
  : Gene ID conversion function using biomart

- [`GetAssayData5()`](https://mengxu98.github.io/scop/reference/GetAssayData5.md)
  :

  Get expression data from `Assay5` or Seurat object

- [`is_outlier()`](https://mengxu98.github.io/scop/reference/is_outlier.md)
  : Detect outliers using MAD (Median Absolute Deviation)

- [`RecoverCounts()`](https://mengxu98.github.io/scop/reference/RecoverCounts.md)
  : Attempt to recover raw counts from the normalized matrix

- [`RenameClusters()`](https://mengxu98.github.io/scop/reference/RenameClusters.md)
  : Rename clusters for the Seurat object

- [`srt_append()`](https://mengxu98.github.io/scop/reference/srt_append.md)
  : Append a Seurat object to another

- [`srt_to_adata()`](https://mengxu98.github.io/scop/reference/srt_to_adata.md)
  : Convert a Seurat object to an AnnData object

- [`srt_reorder()`](https://mengxu98.github.io/scop/reference/srt_reorder.md)
  : Reorder idents by the gene expression

### Feature Processing

- [`AddFeaturesData()`](https://mengxu98.github.io/scop/reference/AddFeaturesData.md)
  : Add features data
- [`AnnotateFeatures()`](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md)
  : Annotate Features
- [`GetFeaturesData()`](https://mengxu98.github.io/scop/reference/GetFeaturesData.md)
  : Get features data
- [`GetSimilarFeatures()`](https://mengxu98.github.io/scop/reference/GetSimilarFeatures.md)
  : Find genes with expression patterns similar to the genes you've
  specified.
- [`RenameFeatures()`](https://mengxu98.github.io/scop/reference/RenameFeatures.md)
  : Rename features for the Seurat object

### Quality Control

- [`db_DoubletDetection()`](https://mengxu98.github.io/scop/reference/db_DoubletDetection.md)
  : Run doublet-calling with DoubletDetection
- [`db_scDblFinder()`](https://mengxu98.github.io/scop/reference/db_scDblFinder.md)
  : Run doublet-calling with scDblFinder
- [`db_scds()`](https://mengxu98.github.io/scop/reference/db_scds.md) :
  Run doublet-calling with scds
- [`db_Scrublet()`](https://mengxu98.github.io/scop/reference/db_Scrublet.md)
  : Run doublet-calling with Scrublet
- [`RunCellQC()`](https://mengxu98.github.io/scop/reference/RunCellQC.md)
  : Run cell-level quality control for single cell RNA-seq data.
- [`RunDoubletCalling()`](https://mengxu98.github.io/scop/reference/RunDoubletCalling.md)
  : Run doublet-calling for single cell RNA-seq data.

### Dimensionality Reduction and Clustering

- [`RunDimReduction()`](https://mengxu98.github.io/scop/reference/RunDimReduction.md)
  : Run dimensionality reduction
- [`RunDM()`](https://mengxu98.github.io/scop/reference/RunDM.md) : Run
  DM (diffusion map)
- [`RunFR()`](https://mengxu98.github.io/scop/reference/RunFR.md) : Run
  Force-Directed Layout (Fruchterman-Reingold algorithm)
- [`RunGLMPCA()`](https://mengxu98.github.io/scop/reference/RunGLMPCA.md)
  : Run GLMPCA (generalized version of principal components analysis)
- [`RunMDS()`](https://mengxu98.github.io/scop/reference/RunMDS.md) :
  Run MDS (multi-dimensional scaling)
- [`RunNMF()`](https://mengxu98.github.io/scop/reference/RunNMF.md) :
  Run NMF (non-negative matrix factorization)
- [`RunUMAP2()`](https://mengxu98.github.io/scop/reference/RunUMAP2.md)
  : Run UMAP (Uniform Manifold Approximation and Projection)
- [`RunPaCMAP()`](https://mengxu98.github.io/scop/reference/RunPaCMAP.md)
  : Run PaCMAP (Pairwise Controlled Manifold Approximation)
- [`RunTriMap()`](https://mengxu98.github.io/scop/reference/RunTriMap.md)
  : Run TriMap (Large-scale Dimensionality Reduction Using Triplets)
- [`RunPHATE()`](https://mengxu98.github.io/scop/reference/RunPHATE.md)
  : Run PHATE (Potential of Heat-diffusion for Affinity-based Trajectory
  Embedding)
- [`RunLargeVis()`](https://mengxu98.github.io/scop/reference/RunLargeVis.md)
  : Run LargeVis (Dimensionality Reduction with a LargeVis-like method)
- [`RunKNNPredict()`](https://mengxu98.github.io/scop/reference/RunKNNPredict.md)
  : Run KNN prediction

### Data Integration and Batch Effect Removal

- [`RunCSSMap()`](https://mengxu98.github.io/scop/reference/RunCSSMap.md)
  : Single-cell reference mapping with CSS method
- [`RunHarmony2()`](https://mengxu98.github.io/scop/reference/RunHarmony2.md)
  : Run Harmony algorithm
- [`RunKNNMap()`](https://mengxu98.github.io/scop/reference/RunKNNMap.md)
  : Single-cell reference mapping with KNN method
- [`RunPCAMap()`](https://mengxu98.github.io/scop/reference/RunPCAMap.md)
  : Single-cell reference mapping with PCA method
- [`RunScmap()`](https://mengxu98.github.io/scop/reference/RunScmap.md)
  : Annotate single cells using scmap.
- [`RunSeuratMap()`](https://mengxu98.github.io/scop/reference/RunSeuratMap.md)
  : Single-cell reference mapping with Seurat method
- [`RunSingleR()`](https://mengxu98.github.io/scop/reference/RunSingleR.md)
  : Annotate single cells using SingleR
- [`RunSymphonyMap()`](https://mengxu98.github.io/scop/reference/RunSymphonyMap.md)
  : Single-cell reference mapping with Symphony method
- [`BBKNN_integrate()`](https://mengxu98.github.io/scop/reference/BBKNN_integrate.md)
  : The BBKNN integration function
- [`CSS_integrate()`](https://mengxu98.github.io/scop/reference/CSS_integrate.md)
  : The CSS integration function
- [`ComBat_integrate()`](https://mengxu98.github.io/scop/reference/ComBat_integrate.md)
  : The ComBat integration function
- [`Conos_integrate()`](https://mengxu98.github.io/scop/reference/Conos_integrate.md)
  : The Conos integration function
- [`Harmony_integrate()`](https://mengxu98.github.io/scop/reference/Harmony_integrate.md)
  : The Harmony integration function
- [`LIGER_integrate()`](https://mengxu98.github.io/scop/reference/LIGER_integrate.md)
  : The LIGER integration function
- [`MNN_integrate()`](https://mengxu98.github.io/scop/reference/MNN_integrate.md)
  : The MNN integration function
- [`Scanorama_integrate()`](https://mengxu98.github.io/scop/reference/Scanorama_integrate.md)
  : The Scanorama integration function
- [`Seurat_integrate()`](https://mengxu98.github.io/scop/reference/Seurat_integrate.md)
  : The Seurat integration function
- [`scVI_integrate()`](https://mengxu98.github.io/scop/reference/scVI_integrate.md)
  : The scVI integration function
- [`Uncorrected_integrate()`](https://mengxu98.github.io/scop/reference/Uncorrected_integrate.md)
  : The Uncorrected integration function
- [`fastMNN_integrate()`](https://mengxu98.github.io/scop/reference/fastMNN_integrate.md)
  : The fastMNN integration function

### Differential Analysis and Enrichment

- [`FindExpressedMarkers()`](https://mengxu98.github.io/scop/reference/FindExpressedMarkers.md)
  : Find Expressed Markers
- [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md)
  : Differential gene test
- [`RunGSEA()`](https://mengxu98.github.io/scop/reference/RunGSEA.md) :
  Perform the enrichment analysis (GSEA) on the genes
- [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md)
  : Perform the enrichment analysis (over-representation) on the genes
- [`RunDynamicEnrichment()`](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)
  : RunDynamicEnrichment
- [`RunProportionTest()`](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
  : Proportion Test

### Trajectory Analysis

- [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md) :
  Run PAGA analysis
- [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md)
  : Run Palantir analysis
- [`RunSlingshot()`](https://mengxu98.github.io/scop/reference/RunSlingshot.md)
  : RunSlingshot
- [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md)
  : Run scVelo workflow
- [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md) :
  Run WOT analysis
- [`RunDynamicFeatures()`](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md)
  : RunDynamicFeatures

### Cell Communication Analysis

- [`RunCellChat()`](https://mengxu98.github.io/scop/reference/RunCellChat.md)
  : Run CellChat

### Visualization Functions

- [`CellChatPlot()`](https://mengxu98.github.io/scop/reference/CellChatPlot.md)
  : Plot CellChat analysis results
- [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
  : Cell Dimensional Plot
- [`CellDimPlot3D()`](https://mengxu98.github.io/scop/reference/CellDimPlot3D.md)
  : 3D-Dimensional reduction plot for cell classification visualization.
- [`CellDensityPlot()`](https://mengxu98.github.io/scop/reference/CellDensityPlot.md)
  : Cell density plot
- [`CellCorHeatmap()`](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md)
  : The Cell Correlation Heatmap
- [`CellStatPlot()`](https://mengxu98.github.io/scop/reference/CellStatPlot.md)
  : Statistical plot of cells
- [`DynamicPlot()`](https://mengxu98.github.io/scop/reference/DynamicPlot.md)
  : Plot dynamic features across pseudotime
- [`DynamicHeatmap()`](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md)
  : Heatmap plot for dynamic features along lineages
- [`EnrichmentPlot()`](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md)
  : Enrichment Plot
- [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
  : Visualize feature values on a 2-dimensional reduction plot
- [`FeatureDimPlot3D()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot3D.md)
  : 3D-Dimensional reduction plot for gene expression visualization.
- [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md)
  : Feature Heatmap
- [`FeatureCorPlot()`](https://mengxu98.github.io/scop/reference/FeatureCorPlot.md)
  : Features correlation plot
- [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)
  : Statistical plot of features
- [`GraphPlot()`](https://mengxu98.github.io/scop/reference/GraphPlot.md)
  : Graph Plot
- [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)
  : The Group Heatmap
- [`GSEAPlot()`](https://mengxu98.github.io/scop/reference/GSEAPlot.md)
  : GSEA Plot
- [`LineagePlot()`](https://mengxu98.github.io/scop/reference/LineagePlot.md)
  : Lineage Plot
- [`PAGAPlot()`](https://mengxu98.github.io/scop/reference/PAGAPlot.md)
  : PAGA plot
- [`ProjectionPlot()`](https://mengxu98.github.io/scop/reference/ProjectionPlot.md)
  : Projection Plot
- [`ProportionTestPlot()`](https://mengxu98.github.io/scop/reference/ProportionTestPlot.md)
  : Proportion Test Plot
- [`StatPlot()`](https://mengxu98.github.io/scop/reference/StatPlot.md)
  : Statistic Plot
- [`TACSPlot()`](https://mengxu98.github.io/scop/reference/TACSPlot.md)
  : Transcript-averaged cell scoring (TACS)
- [`VelocityPlot()`](https://mengxu98.github.io/scop/reference/VelocityPlot.md)
  : Velocity Plot
- [`VolcanoPlot()`](https://mengxu98.github.io/scop/reference/VolcanoPlot.md)
  : Volcano Plot

### Plotting Functions

- [`compute_velocity_on_grid()`](https://mengxu98.github.io/scop/reference/compute_velocity_on_grid.md)
  : Compute velocity on grid
- [`cluster_within_group2()`](https://mengxu98.github.io/scop/reference/cluster_within_group2.md)
  : Cluster within group
- [`print(`*`<scop_logo>`*`)`](https://mengxu98.github.io/scop/reference/print.scop_logo.md)
  : print scop logo

### SCExplorer

- [`CreateDataFile()`](https://mengxu98.github.io/scop/reference/CreateDataFile.md)
  : Create data file
- [`CreateMetaFile()`](https://mengxu98.github.io/scop/reference/CreateMetaFile.md)
  : Create Meta File in HDF5 format from Seurat object
- [`FetchH5()`](https://mengxu98.github.io/scop/reference/FetchH5.md) :
  Fetch data from the hdf5 file
- [`PrepareSCExplorer()`](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md)
  : Prepare Seurat objects for the SCExplorer
- [`RunSCExplorer()`](https://mengxu98.github.io/scop/reference/RunSCExplorer.md)
  : Run SCExplorer

### Database Operations

- [`ListDB()`](https://mengxu98.github.io/scop/reference/ListDB.md) :
  List cached databases
- [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md)
  : Prepare the gene annotation databases

### Data

- [`ifnb_sub`](https://mengxu98.github.io/scop/reference/ifnb_sub.md) :
  A subsetted version of 'ifnb' datasets
- [`panc8_sub`](https://mengxu98.github.io/scop/reference/panc8_sub.md)
  : A subsetted version of human 'panc8' datasets
- [`pancreas_sub`](https://mengxu98.github.io/scop/reference/pancreas_sub.md)
  : A subsetted version of mouse 'pancreas' datasets
- [`ref_scHCL`](https://mengxu98.github.io/scop/reference/ref_scHCL.md)
  [`ref_scMCA`](https://mengxu98.github.io/scop/reference/ref_scHCL.md)
  [`ref_scZCL`](https://mengxu98.github.io/scop/reference/ref_scHCL.md)
  : Reference datasets for cell type annotation in single-cell RNA data
- [`words_excluded`](https://mengxu98.github.io/scop/reference/words_excluded.md)
  : Excluded words in keyword enrichment analysis and extraction
