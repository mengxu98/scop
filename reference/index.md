# Package index

### Spatial and Cellular Omics Analysis Pipeline

- [`scop`](https://mengxu98.github.io/scop/reference/scop-package.md)
  [`scop-package`](https://mengxu98.github.io/scop/reference/scop-package.md)
  : Spatial and single-cell omics analysis pipeline
- [`scop_logo()`](https://mengxu98.github.io/scop/reference/scop_logo.md)
  : scop logo
- [`print(`*`<scop_logo>`*`)`](https://mengxu98.github.io/scop/reference/print.scop_logo.md)
  : print scop logo

### Package and Environment Management

- [`check_python()`](https://mengxu98.github.io/scop/reference/check_python.md)
  : Check and install python packages
- [`env_info()`](https://mengxu98.github.io/scop/reference/env_info.md)
  : Print environment information
- [`env_requirements()`](https://mengxu98.github.io/scop/reference/env_requirements.md)
  : Python environment requirements
- [`ListEnv()`](https://mengxu98.github.io/scop/reference/ListEnv.md) :
  List conda-compatible Python environments
- [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
  : Prepare the python environment
- [`remove_python()`](https://mengxu98.github.io/scop/reference/remove_python.md)
  : Remove Python packages from a conda-compatible Python environment
- [`RemoveEnv()`](https://mengxu98.github.io/scop/reference/RemoveEnv.md)
  : Remove a conda-compatible Python environment

### Feature Processing

- [`AddFeaturesData()`](https://mengxu98.github.io/scop/reference/AddFeaturesData.md)
  : Add features data
- [`AnnotateFeatures()`](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md)
  : Annotate Features
- [`GetFeaturesData()`](https://mengxu98.github.io/scop/reference/GetFeaturesData.md)
  : Get features data
- [`GetSimilarFeatures()`](https://mengxu98.github.io/scop/reference/GetSimilarFeatures.md)
  : Find features with expression patterns similar to provided features
- [`RenameFeatures()`](https://mengxu98.github.io/scop/reference/RenameFeatures.md)
  : Rename features for the Seurat object

### Feature and Statistical Plots

- [`CellCorHeatmap()`](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md)
  : The Cell Correlation Heatmap
- [`CellStatPlot()`](https://mengxu98.github.io/scop/reference/CellStatPlot.md)
  : Statistical plot of cells
- [`FeatureCorPlot()`](https://mengxu98.github.io/scop/reference/FeatureCorPlot.md)
  : Features correlation plot
- [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md)
  : Feature Heatmap
- [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)
  : Statistical plot of features
- [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)
  : The Group Heatmap

### Quality Control

- [`db_DoubletDetection()`](https://mengxu98.github.io/scop/reference/db_DoubletDetection.md)
  : Run doublet-calling with DoubletDetection
- [`db_scDblFinder()`](https://mengxu98.github.io/scop/reference/db_scDblFinder.md)
  : Run doublet-calling with scDblFinder
- [`db_scds()`](https://mengxu98.github.io/scop/reference/db_scds.md) :
  Run doublet-calling with scds
- [`db_Scrublet()`](https://mengxu98.github.io/scop/reference/db_Scrublet.md)
  : Run doublet-calling with Scrublet
- [`RunATACQC()`](https://mengxu98.github.io/scop/reference/RunATACQC.md)
  : Run scATAC quality control metrics
- [`RunCellQC()`](https://mengxu98.github.io/scop/reference/RunCellQC.md)
  : Run cell-level quality control
- [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md)
  : Run ambient RNA decontamination with decontX
- [`RunDoubletCalling()`](https://mengxu98.github.io/scop/reference/RunDoubletCalling.md)
  : Run doublet-calling for single cell RNA-seq data.
- [`RunSpotQC()`](https://mengxu98.github.io/scop/reference/RunSpotQC.md)
  : Run spot-level quality control

### Cell Cycle Analysis

- [`CycGenePrefetch()`](https://mengxu98.github.io/scop/reference/CycGenePrefetch.md)
  : Prefetch cell cycle genes
- [`RunCellCycle()`](https://mengxu98.github.io/scop/reference/RunCellCycle.md)
  : Run cell cycle scoring

### scop Pipeline

- [`standard_scop()`](https://mengxu98.github.io/scop/reference/standard_scop.md)
  : Standard workflow for scop
- [`integration_scop()`](https://mengxu98.github.io/scop/reference/integration_scop.md)
  : The integration workflow

### Metacell Analysis

- [`RunMetaCell()`](https://mengxu98.github.io/scop/reference/RunMetaCell.md)
  : Run metacell partitioning for single-cell data
- [`RunmcRigor()`](https://mengxu98.github.io/scop/reference/RunmcRigor.md)
  : Run mcRigor metacell partition assessment

### Metacell Plots

- [`MetaCellPlot()`](https://mengxu98.github.io/scop/reference/MetaCellPlot.md)
  : Visualize metacell partitions on a dimensionality reduction

### Data Integration and Batch Effect Removal

- [`Uncorrected_integrate()`](https://mengxu98.github.io/scop/reference/Uncorrected_integrate.md)
  : The Uncorrected integration function
- [`Seurat_integrate()`](https://mengxu98.github.io/scop/reference/Seurat_integrate.md)
  : The Seurat integration function
- [`CCA_integrate()`](https://mengxu98.github.io/scop/reference/CCA_integrate.md)
  : Seurat v5 CCA integration
- [`RPCA_integrate()`](https://mengxu98.github.io/scop/reference/RPCA_integrate.md)
  : Seurat v5 RPCA integration
- [`fastMNN_integrate()`](https://mengxu98.github.io/scop/reference/fastMNN_integrate.md)
  : The fastMNN integration function
- [`fastMNN5_integrate()`](https://mengxu98.github.io/scop/reference/fastMNN5_integrate.md)
  : Seurat v5 fastMNN integration
- [`Harmony_integrate()`](https://mengxu98.github.io/scop/reference/Harmony_integrate.md)
  : The Harmony integration function
- [`Harmony5_integrate()`](https://mengxu98.github.io/scop/reference/Harmony5_integrate.md)
  : Seurat v5 Harmony integration
- [`MNN_integrate()`](https://mengxu98.github.io/scop/reference/MNN_integrate.md)
  : The MNN integration function
- [`Scanorama_integrate()`](https://mengxu98.github.io/scop/reference/Scanorama_integrate.md)
  : The Scanorama integration function
- [`BBKNN_integrate()`](https://mengxu98.github.io/scop/reference/BBKNN_integrate.md)
  : The BBKNN integration function
- [`CSS_integrate()`](https://mengxu98.github.io/scop/reference/CSS_integrate.md)
  : The CSS integration function
- [`GLUE_integrate()`](https://mengxu98.github.io/scop/reference/GLUE_integrate.md)
  : The GLUE integration function
- [`LIGER_integrate()`](https://mengxu98.github.io/scop/reference/LIGER_integrate.md)
  : The LIGER integration function
- [`scVI_integrate()`](https://mengxu98.github.io/scop/reference/scVI_integrate.md)
  : The scVI integration function
- [`scVI5_integrate()`](https://mengxu98.github.io/scop/reference/scVI5_integrate.md)
  : Seurat v5 scVI integration
- [`MultiMAP_integrate()`](https://mengxu98.github.io/scop/reference/MultiMAP_integrate.md)
  : The MultiMAP integration function
- [`Conos_integrate()`](https://mengxu98.github.io/scop/reference/Conos_integrate.md)
  : The Conos integration function
- [`ComBat_integrate()`](https://mengxu98.github.io/scop/reference/ComBat_integrate.md)
  : The ComBat integration function
- [`Coralysis_integrate()`](https://mengxu98.github.io/scop/reference/Coralysis_integrate.md)
  : The Coralysis integration function
- [`WNN_integrate()`](https://mengxu98.github.io/scop/reference/WNN_integrate.md)
  : The WNN integration function

### Integration Mapping

- [`RunCSSMap()`](https://mengxu98.github.io/scop/reference/RunCSSMap.md)
  : Single-cell reference mapping with CSS method
- [`RunHarmony2()`](https://mengxu98.github.io/scop/reference/RunHarmony2.md)
  : Run Harmony algorithm
- [`RunLISI()`](https://mengxu98.github.io/scop/reference/RunLISI.md) :
  Compute LISI scores on a Seurat object
- [`RunPCAMap()`](https://mengxu98.github.io/scop/reference/RunPCAMap.md)
  : Single-cell reference mapping with PCA method
- [`RunSeuratMap()`](https://mengxu98.github.io/scop/reference/RunSeuratMap.md)
  : Single-cell reference mapping with Seurat method
- [`RunSymphonyMap()`](https://mengxu98.github.io/scop/reference/RunSymphonyMap.md)
  : Single-cell reference mapping with Symphony method

### Dimensionality Reduction and Clustering

- [`RunDimsEstimate()`](https://mengxu98.github.io/scop/reference/RunDimsEstimate.md)
  : Estimate useful dimensions from a reduction
- [`RunDimsReduction()`](https://mengxu98.github.io/scop/reference/RunDimsReduction.md)
  : Run dimension reduction
- [`RunDM()`](https://mengxu98.github.io/scop/reference/RunDM.md) : Run
  diffusion map (DM)
- [`RunFR()`](https://mengxu98.github.io/scop/reference/RunFR.md) : Run
  Force-Directed Layout (Fruchterman-Reingold algorithm)
- [`RunGLMPCA()`](https://mengxu98.github.io/scop/reference/RunGLMPCA.md)
  : Run generalized principal components analysis (GLMPCA)
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

### Dimensionality Reduction Plots

- [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
  : Cell Dimensional Plot
- [`CellDimPlot3D()`](https://mengxu98.github.io/scop/reference/CellDimPlot3D.md)
  : 3D-Dimensional reduction plot for cell classification visualization.
- [`CellDensityPlot()`](https://mengxu98.github.io/scop/reference/CellDensityPlot.md)
  : Cell density plot
- [`ClusterTreePlot()`](https://mengxu98.github.io/scop/reference/ClusterTreePlot.md)
  : Cluster tree plot
- [`DimsEstimatePlot()`](https://mengxu98.github.io/scop/reference/DimsEstimatePlot.md)
  : Dimension estimate diagnostic plot
- [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
  : Visualize feature values on a 2-dimensional reduction plot
- [`FeatureDimPlot3D()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot3D.md)
  : 3D-Dimensional reduction plot for gene expression visualization.
- [`NMFHeatmap()`](https://mengxu98.github.io/scop/reference/NMFHeatmap.md)
  : NMF similarity heatmap
- [`ProjectionPlot()`](https://mengxu98.github.io/scop/reference/ProjectionPlot.md)
  : Projection Plot
- [`TACSPlot()`](https://mengxu98.github.io/scop/reference/TACSPlot.md)
  : Transcript-averaged cell scoring (TACS)

### Spatial Analysis

- [`SeuratToScopGiotto()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`CreateScopGiotto()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`scop_giotto()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`RunGiottoWorkflow()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`GiottoPreprocess()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`GiottoReduce()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`GiottoCluster()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`GiottoSpatialNetwork()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`GiottoSpatialGenes()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`GiottoSpatialModules()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`GiottoCellProximity()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`GiottoHMRF()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  [`AddGiottoToSeurat()`](https://mengxu98.github.io/scop/reference/SeuratToScopGiotto.md)
  : scop-managed Giotto object workflow
- [`RunBayesSpace()`](https://mengxu98.github.io/scop/reference/RunBayesSpace.md)
  : Run BayesSpace spatial clustering
- [`RunBANKSY()`](https://mengxu98.github.io/scop/reference/RunBANKSY.md)
  : Run BANKSY spatial clustering
- [`RunCARD()`](https://mengxu98.github.io/scop/reference/RunCARD.md) :
  Run CARD spatial deconvolution
- [`RunCytoSPACE()`](https://mengxu98.github.io/scop/reference/RunCytoSPACE.md)
  : Run CytoSPACE spatial assignment
- [`RunSTdeconvolve()`](https://mengxu98.github.io/scop/reference/RunSTdeconvolve.md)
  : Run STdeconvolve reference-free spatial deconvolution
- [`RunGiottoCluster()`](https://mengxu98.github.io/scop/reference/RunGiottoCluster.md)
  : Run Giotto nearest-network clustering
- [`RunGiottoCellProximity()`](https://mengxu98.github.io/scop/reference/RunGiottoCellProximity.md)
  : Run Giotto cell proximity enrichment
- [`RunGiottoSpatialGenes()`](https://mengxu98.github.io/scop/reference/RunGiottoSpatialGenes.md)
  : Run Giotto spatial gene detection
- [`RunGiottoSpatialModules()`](https://mengxu98.github.io/scop/reference/RunGiottoSpatialModules.md)
  : Run Giotto spatial co-expression modules
- [`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md) :
  Run RCTD spatial deconvolution
- [`RunCSIDE()`](https://mengxu98.github.io/scop/reference/RunCSIDE.md)
  : Run C-SIDE spatial differential expression
- [`RunSpaNorm()`](https://mengxu98.github.io/scop/reference/RunSpaNorm.md)
  : Run SpaNorm spatial normalization
- [`RunSmoothClust()`](https://mengxu98.github.io/scop/reference/RunSmoothClust.md)
  : Run smoothclust spatial domain clustering
- [`RunMERINGUE()`](https://mengxu98.github.io/scop/reference/RunMERINGUE.md)
  : Run MERINGUE spatial autocorrelation analysis
- [`RunSpotSweeper()`](https://mengxu98.github.io/scop/reference/RunSpotSweeper.md)
  : Run SpotSweeper spatial quality control
- [`RunSPOTlight()`](https://mengxu98.github.io/scop/reference/RunSPOTlight.md)
  : Run SPOTlight spatial deconvolution
- [`RunSemlaSpatialNetwork()`](https://mengxu98.github.io/scop/reference/RunSemla.md)
  [`RunSemlaLocalG()`](https://mengxu98.github.io/scop/reference/RunSemla.md)
  [`RunSemlaRegionNeighbors()`](https://mengxu98.github.io/scop/reference/RunSemla.md)
  [`RunSemlaRadialDistance()`](https://mengxu98.github.io/scop/reference/RunSemla.md)
  : Run semla spatial analysis backends
- [`RunSpatialEcoTyper()`](https://mengxu98.github.io/scop/reference/RunSpatialEcoTyper.md)
  : Run SpatialEcoTyper spatial ecotype analysis
- [`RunSpatialGradientFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialGradientFeatures.md)
  : Run spatial gradient feature screening
- [`RunSpatialIntegration()`](https://mengxu98.github.io/scop/reference/RunSpatialIntegration.md)
  : Run multi-sample spatial integration
- [`RunSpatialNeighborhood()`](https://mengxu98.github.io/scop/reference/RunSpatialNeighborhood.md)
  : Run spatial neighborhood statistics
- [`RunSpatialVariableFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialVariableFeatures.md)
  : Run spatial variable feature detection

### Spatial Visualization

- [`GiottoPlot()`](https://mengxu98.github.io/scop/reference/GiottoPlot.md)
  : Plot Giotto backend results
- [`SpatialEcoTyperCompositionPlot()`](https://mengxu98.github.io/scop/reference/SpatialEcoTyperCompositionPlot.md)
  : SpatialEcoTyper composition plot
- [`SpatialEcoTyperSpatialPlot()`](https://mengxu98.github.io/scop/reference/SpatialEcoTyperSpatialPlot.md)
  : SpatialEcoTyper spatial plot
- [`SpatialGradientPlot()`](https://mengxu98.github.io/scop/reference/SpatialGradientPlot.md)
  : Plot spatial gradient screening results
- [`SpatialIntegrationPlot()`](https://mengxu98.github.io/scop/reference/SpatialIntegrationPlot.md)
  : Plot spatial integration results
- [`SpatialNeighborhoodPlot()`](https://mengxu98.github.io/scop/reference/SpatialNeighborhoodPlot.md)
  : Spatial neighborhood plot
- [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md)
  : Spatial spot plot
- [`SpatialVariableFeaturePlot()`](https://mengxu98.github.io/scop/reference/SpatialVariableFeaturePlot.md)
  : Plot spatial variable feature results
- [`STdeconvolvePlot()`](https://mengxu98.github.io/scop/reference/STdeconvolvePlot.md)
  : Plot STdeconvolve topic proportions

### Cell Type Annotation

- [`CellTypistModels()`](https://mengxu98.github.io/scop/reference/CellTypistModels.md)
  : Get available CellTypist models
- [`RunCellTypist()`](https://mengxu98.github.io/scop/reference/RunCellTypist.md)
  : Run CellTypist cell type annotation
- [`RunCoEmbedding()`](https://mengxu98.github.io/scop/reference/RunCoEmbedding.md)
  : Co-embed reference and query cells
- [`RunKNNMap()`](https://mengxu98.github.io/scop/reference/RunKNNMap.md)
  : Single-cell reference mapping with KNN method
- [`RunKNNPredict()`](https://mengxu98.github.io/scop/reference/RunKNNPredict.md)
  : Run KNN prediction
- [`RunLabelTransfer()`](https://mengxu98.github.io/scop/reference/RunLabelTransfer.md)
  : Transfer reference labels to query cells
- [`RunReferenceMapping()`](https://mengxu98.github.io/scop/reference/RunReferenceMapping.md)
  : Map query cells into a reference space
- [`RunScmap()`](https://mengxu98.github.io/scop/reference/RunScmap.md)
  : Annotate single cells using scmap.
- [`RunSciBet()`](https://mengxu98.github.io/scop/reference/RunSciBet.md)
  : Annotate single cells using native SciBet
- [`RunSingleR()`](https://mengxu98.github.io/scop/reference/RunSingleR.md)
  : Annotate single cells using SingleR
- [`TrainCellTypist()`](https://mengxu98.github.io/scop/reference/TrainCellTypist.md)
  : Train a CellTypist model

### Differential Expression and Perturbation

- [`FindExpressedMarkers()`](https://mengxu98.github.io/scop/reference/FindExpressedMarkers.md)
  : Find Expressed Markers
- [`RunAugur()`](https://mengxu98.github.io/scop/reference/RunAugur.md)
  : Prioritize perturbed cell types using Augur
- [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md)
  : Differential gene test
- [`RunRareQ()`](https://mengxu98.github.io/scop/reference/RunRareQ.md)
  : RareQ rare-cell population detection
- [`RunScissor()`](https://mengxu98.github.io/scop/reference/RunScissor.md)
  : Run Scissor phenotype-associated cell selection
- [`RunscTenifoldKnk()`](https://mengxu98.github.io/scop/reference/RunScTenifoldKnk.md)
  : Run scTenifoldKnk in-silico knockout analysis
- [`RunscTenifoldNet()`](https://mengxu98.github.io/scop/reference/RunscTenifoldNet.md)
  : Run scTenifoldNet network comparison

### Differential Expression and Perturbation Plots

- [`DEtestPlot()`](https://mengxu98.github.io/scop/reference/DEtestPlot.md)
  : Differential Expression Test Plot
- [`DEtestManhattanPlot()`](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md)
  : DEtest Manhattan Plot
- [`DEtestRingPlot()`](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)
  : DEtest Ring Plot
- [`VolcanoPlot()`](https://mengxu98.github.io/scop/reference/VolcanoPlot.md)
  : Volcano Plot
- [`ScissorPlot()`](https://mengxu98.github.io/scop/reference/ScissorPlot.md)
  : Plot Scissor results
- [`scTenifoldKnkPlot()`](https://mengxu98.github.io/scop/reference/scTenifoldKnkPlot.md)
  : scTenifoldKnk Plot
- [`scTenifoldNetPlot()`](https://mengxu98.github.io/scop/reference/scTenifoldNetPlot.md)
  : scTenifoldNet Plot

### Bulk and Composition Analysis

- [`RunDeconvolution()`](https://mengxu98.github.io/scop/reference/RunDeconvolution.md)
  : Run bulk or pseudobulk deconvolution
- [`RunCIBERSORT()`](https://mengxu98.github.io/scop/reference/RunCIBERSORT.md)
  : Run CIBERSORT deconvolution
- [`RunESTIMATE()`](https://mengxu98.github.io/scop/reference/RunESTIMATE.md)
  : Run ESTIMATE tumor microenvironment scoring
- [`RunMilo()`](https://mengxu98.github.io/scop/reference/RunMilo.md) :
  Milo differential abundance wrapper
- [`RunPermutation()`](https://mengxu98.github.io/scop/reference/RunPermutation.md)
  : Permutation-based proportion test
- [`RunPropeller()`](https://mengxu98.github.io/scop/reference/RunPropeller.md)
  : Propeller differential abundance wrapper
- [`RunProportionTest()`](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
  : Proportion Test
- [`RunscCODA()`](https://mengxu98.github.io/scop/reference/RunscCODA.md)
  : scCODA differential abundance

### Bulk and Composition Plots

- [`DeconvolutionPlot()`](https://mengxu98.github.io/scop/reference/DeconvolutionPlot.md)
  : Plot deconvolution results
- [`ImmuneAbundancePlot()`](https://mengxu98.github.io/scop/reference/ImmuneAbundancePlot.md)
  : Immune abundance plots
- [`EstimateScorePlot()`](https://mengxu98.github.io/scop/reference/EstimateScorePlot.md)
  : ESTIMATE score plots
- [`EstimateGenePlot()`](https://mengxu98.github.io/scop/reference/EstimateGenePlot.md)
  : Gene and ESTIMATE score relationship plots
- [`GeneImmuneCorPlot()`](https://mengxu98.github.io/scop/reference/GeneImmuneCorPlot.md)
  : Gene-immune correlation butterfly plot
- [`ProportionTestPlot()`](https://mengxu98.github.io/scop/reference/ProportionTestPlot.md)
  : Proportion Test Plot

### Copy Number Variation Analysis

- [`RunCNV()`](https://mengxu98.github.io/scop/reference/RunCNV.md) :
  Run copy-number alteration inference
- [`CNVPlot()`](https://mengxu98.github.io/scop/reference/CNVPlot.md) :
  Plot SCOP CNV results

### Enrichment Analysis

- [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md)
  : Cell scoring
- [`RunDorothea()`](https://mengxu98.github.io/scop/reference/RunDorothea.md)
  : Run DoRothEA transcription factor activity inference
- [`RunDynamicEnrichment()`](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)
  : RunDynamicEnrichment
- [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md)
  : Perform the enrichment analysis (over-representation) on the genes
- [`RunGSEA()`](https://mengxu98.github.io/scop/reference/RunGSEA.md) :
  Perform the enrichment analysis (GSEA) on the genes
- [`RunGSVA()`](https://mengxu98.github.io/scop/reference/RunGSVA.md) :
  Perform Gene Set Variation Analysis (GSVA)
- [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md)
  : Run metabolism pathway scoring

### Enrichment Analysis Plot

- [`DorotheaPlot()`](https://mengxu98.github.io/scop/reference/DorotheaPlot.md)
  : Plot differential DoRothEA TF activity
- [`EnrichmentPlot()`](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md)
  : Enrichment Plot
- [`FerrisWheelPlot()`](https://mengxu98.github.io/scop/reference/FerrisWheelPlot.md)
  : Ferris Wheel Plot
- [`GSEAPlot()`](https://mengxu98.github.io/scop/reference/GSEAPlot.md)
  : GSEA Plot
- [`GSVAPlot()`](https://mengxu98.github.io/scop/reference/GSVAPlot.md)
  : Plots for GSVA (Gene Set Variation Analysis)
- [`MetabolismPlot()`](https://mengxu98.github.io/scop/reference/MetabolismPlot.md)
  : Plots for metabolism pathway scoring

### Metabolic Flux Analysis

- [`RunscFEA()`](https://mengxu98.github.io/scop/reference/RunscFEA.md)
  : Run scFEA flux estimation for a Seurat object
- [`scFEAHeatmap()`](https://mengxu98.github.io/scop/reference/scFEAHeatmap.md)
  : Plot scFEA module flux heatmap
- [`scFEAVolcanoPlot()`](https://mengxu98.github.io/scop/reference/scFEAVolcanoPlot.md)
  : Plot scFEA flux Cohen's d volcano plots
- [`scFEABalanceBarPlot()`](https://mengxu98.github.io/scop/reference/scFEABalanceBarPlot.md)
  : Plot scFEA metabolite balance changes

### Dynamic Trajectory Analysis

- [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md)
  : Run CellRank analysis
- [`RunCytoTRACE()`](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md)
  : Run CytoTRACE 2
- [`RunMonocle2()`](https://mengxu98.github.io/scop/reference/RunMonocle2.md)
  : Run Monocle2 analysis
- [`RunMonocle3()`](https://mengxu98.github.io/scop/reference/RunMonocle3.md)
  : Run Monocle3 analysis
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
  : Calculates dynamic features for lineages

### Dynamic Trajectory Plots

- [`BranchStreamPlot()`](https://mengxu98.github.io/scop/reference/BranchStreamPlot.md)
  : Branch Stream Plot
- [`CytoTRACEPlot()`](https://mengxu98.github.io/scop/reference/CytoTRACEPlot.md)
  : Plot CytoTRACE 2 Results
- [`DynamicPlot()`](https://mengxu98.github.io/scop/reference/DynamicPlot.md)
  : Plot dynamic features across pseudotime
- [`DynamicHeatmap()`](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md)
  : Heatmap plot for dynamic features along lineages
- [`LineagePlot()`](https://mengxu98.github.io/scop/reference/LineagePlot.md)
  : Lineage Plot
- [`PAGAPlot()`](https://mengxu98.github.io/scop/reference/PAGAPlot.md)
  : PAGA plot
- [`PalantirTrajectoryPlot()`](https://mengxu98.github.io/scop/reference/PalantirTrajectoryPlot.md)
  : Plot Palantir trajectories
- [`PseudotimeProjectionPlot()`](https://mengxu98.github.io/scop/reference/PseudotimeProjectionPlot.md)
  : Pseudotime Projection Plot
- [`VelocityPlot()`](https://mengxu98.github.io/scop/reference/VelocityPlot.md)
  : Velocity Plot

### Cell-Cell Communication Analysis

- [`RunCCC()`](https://mengxu98.github.io/scop/reference/RunCCC.md) :
  Run common cell-cell communication analyses
- [`RunCellChat()`](https://mengxu98.github.io/scop/reference/RunCellChat.md)
  : Run CellChat analysis
- [`RunCellphoneDB()`](https://mengxu98.github.io/scop/reference/RunCellphoneDB.md)
  : Run CellphoneDB analysis
- [`RunLIANA()`](https://mengxu98.github.io/scop/reference/RunLIANA.md)
  : Run LIANA cell-cell communication analysis
- [`RunNichenetr()`](https://mengxu98.github.io/scop/reference/RunNichenetr.md)
  : Run NicheNet analysis
- [`RunMultiNichenetr()`](https://mengxu98.github.io/scop/reference/RunMultiNichenetr.md)
  : Run MultiNicheNet analysis
- [`RunscOMM()`](https://mengxu98.github.io/scop/reference/RunscOMM.md)
  : Run scOMM label prediction
- [`ccc_to_adata()`](https://mengxu98.github.io/scop/reference/ccc_to_adata.md)
  : Convert CCC results to OmicVerse communication AnnData
- [`ccc_to_liana()`](https://mengxu98.github.io/scop/reference/ccc_to_liana.md)
  : Convert CCC results to a LIANA-like table

### Cell-Cell Communication Plots

- [`CCCHeatmap()`](https://mengxu98.github.io/scop/reference/CCCHeatmap.md)
  : CCC heatmap and dot matrix plot
- [`CCCNetworkPlot()`](https://mengxu98.github.io/scop/reference/CCCNetworkPlot.md)
  : CCC network and flow plots
- [`CCCStatPlot()`](https://mengxu98.github.io/scop/reference/CCCStatPlot.md)
  : CCC statistical distribution and summary plots

### Gene Regulatory Network Inference

- [`RunGRNBoost2()`](https://mengxu98.github.io/scop/reference/RunGRNBoost2.md)
  : Infer gene regulatory networks with GRNBoost2
- [`RunGENIE3()`](https://mengxu98.github.io/scop/reference/RunGENIE3.md)
  : Infer gene regulatory networks with GENIE3
- [`RunSCENIC()`](https://mengxu98.github.io/scop/reference/RunSCENIC.md)
  : Run SCENIC gene regulatory network analysis
- [`RunSCENICPlus()`](https://mengxu98.github.io/scop/reference/RunSCENICPlus.md)
  : Run SCENICPlus-style eGRN analysis

### Gene Regulatory Network Plots

- [`SCENICPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlot.md)
  : Plot top regulon specificity scores from SCENIC results

### Other

- [`BenchmarkPlot()`](https://mengxu98.github.io/scop/reference/BenchmarkPlot.md)
  : Plot benchmark metrics
- [`CoverageTrackPlot()`](https://mengxu98.github.io/scop/reference/CoverageTrackPlot.md)
  : Coverage track plot for ATAC data
- [`LISIPlot()`](https://mengxu98.github.io/scop/reference/LISIPlot.md)
  : Plot LISI scores
- [`tAgePlot()`](https://mengxu98.github.io/scop/reference/tAgePlot.md)
  : Plot tAge transcriptomic aging-clock predictions
- [`RuntAge()`](https://mengxu98.github.io/scop/reference/RuntAge.md) :
  Run tAge transcriptomic aging-clock prediction

### SCExplorer

- [`CreateDataFile()`](https://mengxu98.github.io/scop/reference/CreateDataFile.md)
  : Create HDF5 data file from Seurat object
- [`CreateMetaFile()`](https://mengxu98.github.io/scop/reference/CreateMetaFile.md)
  : Create Meta File in HDF5 format from Seurat object
- [`FetchH5()`](https://mengxu98.github.io/scop/reference/FetchH5.md) :
  Fetch data from the hdf5 file and returns a Seurat object
- [`PrepareSCExplorer()`](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md)
  : Prepare Seurat objects for the SCExplorer
- [`RunSCExplorer()`](https://mengxu98.github.io/scop/reference/RunSCExplorer.md)
  : Run SCExplorer

### Data Conversion

- [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md)
  : Convert an anndata object to a seurat object

- [`h5ad_to_srt()`](https://mengxu98.github.io/scop/reference/h5ad_to_srt.md)
  :

  Read an `.h5ad` file and convert to a `Seurat`

- [`loom_to_adata()`](https://mengxu98.github.io/scop/reference/loom_to_adata.md)
  :

  Read a `.loom` file as an AnnData object

- [`loom_to_srt()`](https://mengxu98.github.io/scop/reference/loom_to_srt.md)
  :

  Read a `.loom` file and convert to a `Seurat`

- [`srt_to_adata()`](https://mengxu98.github.io/scop/reference/srt_to_adata.md)
  : Convert a Seurat object to an AnnData object

- [`srt_to_h5ad()`](https://mengxu98.github.io/scop/reference/srt_to_h5ad.md)
  :

  Convert a Seurat object to an `.h5ad` file

### Data Processing

- [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  :

  Check and report the type of data in `Seurat` object

- [`CheckDataList()`](https://mengxu98.github.io/scop/reference/CheckDataList.md)
  :

  Check and preprocess a list of `Seurat` objects

- [`CheckDataMerge()`](https://mengxu98.github.io/scop/reference/CheckDataMerge.md)
  : Check and preprocess a merged seurat object

- [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
  : Find the default reduction name in a Seurat object

- [`FetchDataZero()`](https://mengxu98.github.io/scop/reference/FetchDataZero.md)
  : FetchData but with zeroes for unavailable genes

- [`GeneConvert()`](https://mengxu98.github.io/scop/reference/GeneConvert.md)
  : Gene ID conversion function using biomart

- [`ConvertHomologs()`](https://mengxu98.github.io/scop/reference/ConvertHomologs.md)
  : Convert homologous gene symbols in expression objects

- [`GetAssayData5()`](https://mengxu98.github.io/scop/reference/GetAssayData5.md)
  :

  Get expression data from `Assay5` or Seurat object

- [`RecoverCounts()`](https://mengxu98.github.io/scop/reference/RecoverCounts.md)
  : Attempt to recover raw counts from the normalized matrix

- [`RenameClusters()`](https://mengxu98.github.io/scop/reference/RenameClusters.md)
  : Rename clusters for the Seurat object

- [`srt_append()`](https://mengxu98.github.io/scop/reference/srt_append.md)
  : Append a Seurat object to another

- [`srt_reorder()`](https://mengxu98.github.io/scop/reference/srt_reorder.md)
  : Reorder idents by the gene expression

### Database Operations

- [`ListDB()`](https://mengxu98.github.io/scop/reference/ListDB.md) :
  List cached databases
- [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md)
  : Prepare the gene annotation databases
- [`RunCisTarget()`](https://mengxu98.github.io/scop/reference/RunCisTarget.md)
  : Run cisTarget motif enrichment on a GRN adjacency table

### Data

- [`islet_bulk`](https://mengxu98.github.io/scop/reference/islet_bulk.md)
  : Human pancreatic islet bulk RNA-seq example dataset
- [`panc8_sub`](https://mengxu98.github.io/scop/reference/panc8_sub.md)
  : A subsetted version of human 'panc8' datasets
- [`pancreas_sub`](https://mengxu98.github.io/scop/reference/pancreas_sub.md)
  : A subsetted version of mouse 'pancreas' datasets
- [`pbmcmultiome_sub`](https://mengxu98.github.io/scop/reference/pbmcmultiome_sub.md)
  : A small human PBMC multiome example dataset
- [`ref_scMCA`](https://mengxu98.github.io/scop/reference/ref_scMCA.md)
  : Reference datasets for cell type annotation in single-cell RNA data
- [`visium_human_pancreas_sub`](https://mengxu98.github.io/scop/reference/visium_human_pancreas_sub.md)
  : A human pancreas Visium spatial example dataset
- [`words_excluded`](https://mengxu98.github.io/scop/reference/words_excluded.md)
  : Excluded words in keyword enrichment analysis and extraction
