## Overview and Workflows

### Package and Standard Workflows

- [`scop`](https://mengxu98.github.io/scop/reference/scop-package.md)
  [`scop-package`](https://mengxu98.github.io/scop/reference/scop-package.md)
  : Single-cell and Spatial omics analysis pipeline
- [`scop_logo()`](https://mengxu98.github.io/scop/reference/scop_logo.md)
  [`print(`*`<scop_logo>`*`)`](https://mengxu98.github.io/scop/reference/scop_logo.md)
  : scop logo
- [`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md)
  : Standard single-cell and spatial workflow

## Quality Control and Preprocessing

### Cell Quality and Contamination Control

- [`RunCellQC()`](https://mengxu98.github.io/scop/reference/RunCellQC.md)
  : Cell-level quality control
- [`RunDoubletCalling()`](https://mengxu98.github.io/scop/reference/RunDoubletCalling.md)
  : Doublet calling
- [`RunDoubletDetection()`](https://mengxu98.github.io/scop/reference/RunDoubletDetection.md)
  : Run doublet-calling with DoubletDetection
- [`RunscDblFinder()`](https://mengxu98.github.io/scop/reference/RunscDblFinder.md)
  : Run doublet-calling with scDblFinder
- [`Runscds()`](https://mengxu98.github.io/scop/reference/Runscds.md) :
  Run doublet-calling with scds
- [`RunScrublet()`](https://mengxu98.github.io/scop/reference/RunScrublet.md)
  : Run doublet-calling with Scrublet
- [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md)
  : Ambient RNA decontamination with decontX

### Normalization and Feature Selection

- [`NormalizeData()`](https://mengxu98.github.io/scop/reference/NormalizeData.md)
  : Normalize a single-cell object
- [`FindVariableFeatures()`](https://mengxu98.github.io/scop/reference/FindVariableFeatures.md)
  : Find variable features
- [`ScaleData()`](https://mengxu98.github.io/scop/reference/ScaleData.md)
  : Scale expression data

### Cell Cycle Analysis

- [`CellCycleScoring()`](https://mengxu98.github.io/scop/reference/CellCycleScoring.md)
  : Cell-cycle scoring with native sparse module scoring
- [`CycGenePrefetch()`](https://mengxu98.github.io/scop/reference/CycGenePrefetch.md)
  : Prefetch cell cycle genes
- [`RunCellCycle()`](https://mengxu98.github.io/scop/reference/RunCellCycle.md)
  : Run cell cycle scoring

## Dimensionality Reduction and Clustering

### Reduction Selection and Linear Embeddings

- [`RunDimsEstimate()`](https://mengxu98.github.io/scop/reference/RunDimsEstimate.md)
  : Estimate useful dimensions from a reduction
- [`DimsEstimatePlot()`](https://mengxu98.github.io/scop/reference/DimsEstimatePlot.md)
  : Dimension estimate diagnostic plot
- [`RunDimsReduction()`](https://mengxu98.github.io/scop/reference/RunDimsReduction.md)
  : Run dimension reduction
- [`RunPCA()`](https://mengxu98.github.io/scop/reference/RunPCA.md) :
  Run principal component analysis
- [`RunCCA()`](https://mengxu98.github.io/scop/reference/RunCCA.md) :
  Run canonical correlation analysis
- [`RunGLMPCA()`](https://mengxu98.github.io/scop/reference/RunGLMPCA.md)
  : Run generalized principal components analysis (GLMPCA)
- [`RunMDS()`](https://mengxu98.github.io/scop/reference/RunMDS.md) :
  Run MDS (multi-dimensional scaling)
- [`RunNMF()`](https://mengxu98.github.io/scop/reference/RunNMF.md) :
  Run NMF
- [`NMFHeatmap()`](https://mengxu98.github.io/scop/reference/NMFHeatmap.md)
  : NMF similarity heatmap

### Nonlinear Embeddings

- [`RunUMAP()`](https://mengxu98.github.io/scop/reference/RunUMAP.md) :
  Run UMAP
- [`RunUMAP2()`](https://mengxu98.github.io/scop/reference/RunUMAP2.md)
  : Run UMAP
- [`RunDM()`](https://mengxu98.github.io/scop/reference/RunDM.md) : Run
  diffusion map (DM)
- [`RunFR()`](https://mengxu98.github.io/scop/reference/RunFR.md) : Run
  Force-Directed Layout (Fruchterman-Reingold algorithm)
- [`RunPaCMAP()`](https://mengxu98.github.io/scop/reference/RunPaCMAP.md)
  : Run PaCMAP
- [`RunTriMap()`](https://mengxu98.github.io/scop/reference/RunTriMap.md)
  : Run TriMap
- [`RunPHATE()`](https://mengxu98.github.io/scop/reference/RunPHATE.md)
  : Run PHATE
- [`RunLargeVis()`](https://mengxu98.github.io/scop/reference/RunLargeVis.md)
  : Run LargeVis (Dimensionality Reduction with a LargeVis-like method)

### Clustering and Neighborhoods

- [`FindNeighbors()`](https://mengxu98.github.io/scop/reference/FindNeighbors.md)
  : Find nearest neighbors
- [`FindClusters()`](https://mengxu98.github.io/scop/reference/FindClusters.md)
  : Seurat clustering with transparent backend delegation
- [`ClusterTreePlot()`](https://mengxu98.github.io/scop/reference/ClusterTreePlot.md)
  : Cluster tree plot
- [`RunCHOIR()`](https://mengxu98.github.io/scop/reference/RunCHOIR.md)
  : Run CHOIR clustering
- [`RunRareQ()`](https://mengxu98.github.io/scop/reference/RunRareQ.md)
  : RareQ rare-cell population detection

### Embedding Visualization

- [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
  : Cell Dimensional Plot
- [`CellDimPlot3D()`](https://mengxu98.github.io/scop/reference/CellDimPlot3D.md)
  : 3D-Dimensional reduction plot for cell classification visualization.
- [`CellDensityPlot()`](https://mengxu98.github.io/scop/reference/CellDensityPlot.md)
  : Cell density plot
- [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
  : Feature values on a 2D reduction
- [`FeatureDimPlot3D()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot3D.md)
  : Feature values on a 3D reduction
- [`ProjectionPlot()`](https://mengxu98.github.io/scop/reference/ProjectionPlot.md)
  : Projection Plot
- [`TACSPlot()`](https://mengxu98.github.io/scop/reference/TACSPlot.md)
  : Transcript-averaged cell scoring (TACS)

## Data Integration and Reference Mapping

### Integration Workflows and Feature Preparation

- [`RunIntegration()`](https://mengxu98.github.io/scop/reference/RunIntegration.md)
  : Integration workflow
- [`SelectIntegrationFeatures()`](https://mengxu98.github.io/scop/reference/SelectIntegrationFeatures.md)
  : Select integration features with exact vectorized ranking
- [`Uncorrected_integrate()`](https://mengxu98.github.io/scop/reference/Uncorrected_integrate.md)
  : Uncorrected integration function

### Anchor and Neighbor-Based Integration

- [`Seurat_integrate()`](https://mengxu98.github.io/scop/reference/Seurat_integrate.md)
  : Seurat integration function
- [`CCA_integrate()`](https://mengxu98.github.io/scop/reference/CCA_integrate.md)
  : Seurat v5 CCA integration
- [`RPCA_integrate()`](https://mengxu98.github.io/scop/reference/RPCA_integrate.md)
  : Seurat v5 RPCA integration
- [`MNN_integrate()`](https://mengxu98.github.io/scop/reference/MNN_integrate.md)
  : MNN integration function
- [`fastMNN_integrate()`](https://mengxu98.github.io/scop/reference/fastMNN_integrate.md)
  : FastMNN integration function
- [`fastMNN5_integrate()`](https://mengxu98.github.io/scop/reference/fastMNN5_integrate.md)
  : Seurat v5 fastMNN integration
- [`Scanorama_integrate()`](https://mengxu98.github.io/scop/reference/Scanorama_integrate.md)
  : Scanorama integration function

### Embedding and Factor Integration

- [`Harmony_integrate()`](https://mengxu98.github.io/scop/reference/Harmony_integrate.md)
  : Harmony integration function
- [`Harmony5_integrate()`](https://mengxu98.github.io/scop/reference/Harmony5_integrate.md)
  : Seurat v5 Harmony integration
- [`RunHarmony2()`](https://mengxu98.github.io/scop/reference/RunHarmony2.md)
  : Run Harmony algorithm
- [`LIGER_integrate()`](https://mengxu98.github.io/scop/reference/LIGER_integrate.md)
  : LIGER integration function
- [`scVI_integrate()`](https://mengxu98.github.io/scop/reference/scVI_integrate.md)
  : ScVI integration function
- [`scVI5_integrate()`](https://mengxu98.github.io/scop/reference/scVI5_integrate.md)
  : Seurat v5 scVI integration
- [`Coralysis_integrate()`](https://mengxu98.github.io/scop/reference/Coralysis_integrate.md)
  : Coralysis integration function

### Graph and Expression-Space Integration

- [`BBKNN_integrate()`](https://mengxu98.github.io/scop/reference/BBKNN_integrate.md)
  : BBKNN integration function
- [`Conos_integrate()`](https://mengxu98.github.io/scop/reference/Conos_integrate.md)
  : Conos integration function
- [`CSS_integrate()`](https://mengxu98.github.io/scop/reference/CSS_integrate.md)
  : CSS integration function
- [`ComBat_integrate()`](https://mengxu98.github.io/scop/reference/ComBat_integrate.md)
  : ComBat integration function

### Multimodal Integration

- [`WNN_integrate()`](https://mengxu98.github.io/scop/reference/WNN_integrate.md)
  : WNN integration function
- [`FindMultiModalNeighbors()`](https://mengxu98.github.io/scop/reference/FindMultiModalNeighbors.md)
  : Seurat weighted-nearest-neighbor construction
- [`GLUE_integrate()`](https://mengxu98.github.io/scop/reference/GLUE_integrate.md)
  : GLUE integration function
- [`MultiMAP_integrate()`](https://mengxu98.github.io/scop/reference/MultiMAP_integrate.md)
  : MultiMAP integration function

### Reference Mapping and Label Transfer

- [`RunReferenceMapping()`](https://mengxu98.github.io/scop/reference/RunReferenceMapping.md)
  : Map query cells into a reference space
- [`RunLabelTransfer()`](https://mengxu98.github.io/scop/reference/RunLabelTransfer.md)
  : Transfer reference labels to query cells
- [`RunSeuratMap()`](https://mengxu98.github.io/scop/reference/RunSeuratMap.md)
  : Single-cell reference mapping with Seurat method
- [`RunSymphonyMap()`](https://mengxu98.github.io/scop/reference/RunSymphonyMap.md)
  : Single-cell reference mapping with Symphony method
- [`RunPCAMap()`](https://mengxu98.github.io/scop/reference/RunPCAMap.md)
  : Single-cell reference mapping with PCA method
- [`RunCSSMap()`](https://mengxu98.github.io/scop/reference/RunCSSMap.md)
  : Single-cell reference mapping with CSS method
- [`RunKNNMap()`](https://mengxu98.github.io/scop/reference/RunKNNMap.md)
  : Single-cell reference mapping with KNN method
- [`RunCoEmbedding()`](https://mengxu98.github.io/scop/reference/RunCoEmbedding.md)
  : Co-embed reference and query cells

### Integration Evaluation

- [`RunIntegrationBenchmark()`](https://mengxu98.github.io/scop/reference/RunIntegrationBenchmark.md)
  : Benchmark single-cell integration methods
- [`IntegrationBenchmarkPlot()`](https://mengxu98.github.io/scop/reference/IntegrationBenchmarkPlot.md)
  : Plot integration benchmark results
- [`RunLISI()`](https://mengxu98.github.io/scop/reference/RunLISI.md) :
  Compute LISI scores on a Seurat object

## Cell Type Annotation

### Reference-Based Annotation and Models

- [`CellTypistModels()`](https://mengxu98.github.io/scop/reference/CellTypistModels.md)
  : Get available CellTypist models
- [`TrainCellTypist()`](https://mengxu98.github.io/scop/reference/TrainCellTypist.md)
  : Train a CellTypist model
- [`RunCellTypist()`](https://mengxu98.github.io/scop/reference/RunCellTypist.md)
  : Run CellTypist cell type annotation
- [`RunSingleR()`](https://mengxu98.github.io/scop/reference/RunSingleR.md)
  : Annotate single cells using SingleR
- [`RunSciBet()`](https://mengxu98.github.io/scop/reference/RunSciBet.md)
  : Annotate single cells with a C++ SciBet-compatible classifier
- [`RunScmap()`](https://mengxu98.github.io/scop/reference/RunScmap.md)
  : Annotate single cells using scmap
- [`RunKNNPredict()`](https://mengxu98.github.io/scop/reference/RunKNNPredict.md)
  : Run KNN prediction

### Specialized Cell State Annotation

- [`RunscMalignantFinder()`](https://mengxu98.github.io/scop/reference/RunscMalignantFinder.md)
  : Run scMalignantFinder malignant cell identification
- [`RunscMalignantRegion()`](https://mengxu98.github.io/scop/reference/RunscMalignantRegion.md)
  : Run scMalignantFinder malignant spatial region identification
- [`RunscMalignantStates()`](https://mengxu98.github.io/scop/reference/RunscMalignantStates.md)
  : Run scMalignantFinder cancer cell state scoring

## Differential Expression and Cell Composition

### Markers and Differential Expression

- [`FindAllMarkers()`](https://mengxu98.github.io/scop/reference/FindAllMarkers.md)
  : Find markers for all groups
- [`FindMarkers()`](https://mengxu98.github.io/scop/reference/FindMarkers.md)
  : Find markers between groups
- [`FindExpressedMarkers()`](https://mengxu98.github.io/scop/reference/FindExpressedMarkers.md)
  : Find Expressed Markers
- [`FoldChange()`](https://mengxu98.github.io/scop/reference/FoldChange.md)
  : Calculate expression fold change
- [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md)
  : Differential gene test
- [`DEtestPlot()`](https://mengxu98.github.io/scop/reference/DEtestPlot.md)
  : Differential Expression Test Plot
- [`DEtestManhattanPlot()`](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md)
  : DEtest Manhattan Plot
- [`DEtestRingPlot()`](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)
  : DEtest Ring Plot
- [`VolcanoPlot()`](https://mengxu98.github.io/scop/reference/VolcanoPlot.md)
  : Volcano Plot

### Expression Aggregation and Metacells

- [`AggregateExpression()`](https://mengxu98.github.io/scop/reference/AggregateExpression.md)
  : Aggregate expression by cell group
- [`AverageExpression()`](https://mengxu98.github.io/scop/reference/AverageExpression.md)
  : Average expression by cell group
- [`RunMetaCell()`](https://mengxu98.github.io/scop/reference/RunMetaCell.md)
  : Run metacell partitioning for single-cell data
- [`MetaCellPlot()`](https://mengxu98.github.io/scop/reference/MetaCellPlot.md)
  : Visualize metacell partitions on a dimensionality reduction
- [`RunmcRigor()`](https://mengxu98.github.io/scop/reference/RunmcRigor.md)
  : Run mcRigor metacell partition assessment

### Cell Composition and Abundance

- [`RunMilo()`](https://mengxu98.github.io/scop/reference/RunMilo.md) :
  Milo differential abundance wrapper
- [`RunPropeller()`](https://mengxu98.github.io/scop/reference/RunPropeller.md)
  : Propeller differential abundance wrapper
- [`RunscCODA()`](https://mengxu98.github.io/scop/reference/RunscCODA.md)
  : scCODA differential abundance
- [`RunProportionTest()`](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
  : Proportion Test
- [`ProportionTestPlot()`](https://mengxu98.github.io/scop/reference/ProportionTestPlot.md)
  : Proportion Test Plot
- [`RunPermutation()`](https://mengxu98.github.io/scop/reference/RunPermutation.md)
  : Permutation-based proportion test

### Bulk Deconvolution and Immune Scores

- [`RunDeconvolution()`](https://mengxu98.github.io/scop/reference/RunDeconvolution.md)
  : Run bulk or pseudobulk deconvolution
- [`DeconvolutionPlot()`](https://mengxu98.github.io/scop/reference/DeconvolutionPlot.md)
  : Plot deconvolution results
- [`RunCIBERSORT()`](https://mengxu98.github.io/scop/reference/RunCIBERSORT.md)
  : Run CIBERSORT deconvolution
- [`ImmuneAbundancePlot()`](https://mengxu98.github.io/scop/reference/ImmuneAbundancePlot.md)
  : Immune abundance plots
- [`RunESTIMATE()`](https://mengxu98.github.io/scop/reference/RunESTIMATE.md)
  : Run ESTIMATE tumor microenvironment scoring
- [`EstimateScorePlot()`](https://mengxu98.github.io/scop/reference/EstimateScorePlot.md)
  : ESTIMATE score plots
- [`EstimateGenePlot()`](https://mengxu98.github.io/scop/reference/EstimateGenePlot.md)
  : Gene and ESTIMATE score relationship plots
- [`GeneImmuneCorPlot()`](https://mengxu98.github.io/scop/reference/GeneImmuneCorPlot.md)
  : Gene-immune correlation butterfly plot
- [`estimate_signatures`](https://mengxu98.github.io/scop/reference/estimate_signatures.md)
  : ESTIMATE gene signatures

## Functional and Pathway Analysis

### Enrichment Testing and Visualization

- [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md)
  : Perform the enrichment analysis (over-representation) on the genes
- [`EnrichmentPlot()`](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md)
  : Enrichment Plot
- [`FerrisWheelPlot()`](https://mengxu98.github.io/scop/reference/FerrisWheelPlot.md)
  : Ferris Wheel Plot
- [`RunGSEA()`](https://mengxu98.github.io/scop/reference/RunGSEA.md) :
  Perform the enrichment analysis (GSEA) on the genes
- [`GSEAPlot()`](https://mengxu98.github.io/scop/reference/GSEAPlot.md)
  : GSEA Plot

### Gene Set and Pathway Scoring

- [`AddModuleScore()`](https://mengxu98.github.io/scop/reference/AddModuleScore.md)
  : Add expression-program scores with a sparse native fast path
- [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md)
  : Cell scoring
- [`CellScoringPlot()`](https://mengxu98.github.io/scop/reference/CellScoringPlot.md)
  : Plot CellScoring results from a Seurat object
- [`RunGSVA()`](https://mengxu98.github.io/scop/reference/RunGSVA.md) :
  Perform Gene Set Variation Analysis (GSVA)
- [`GSVAPlot()`](https://mengxu98.github.io/scop/reference/GSVAPlot.md)
  : Plots for GSVA (Gene Set Variation Analysis)

### Immune Response Programs

- [`PrepareIREA()`](https://mengxu98.github.io/scop/reference/PrepareIREA.md)
  : Prepare Immune Dictionary references
- [`RunIREA()`](https://mengxu98.github.io/scop/reference/RunIREA.md) :
  Analyse cytokine responses or immune-cell polarization
- [`IREAPlot()`](https://mengxu98.github.io/scop/reference/IREAPlot.md)
  : Plot IREA response or polarization results

### Metabolic Activity and Flux

- [`RunMetabolism()`](https://mengxu98.github.io/scop/reference/RunMetabolism.md)
  : Run metabolism pathway scoring
- [`MetabolismPlot()`](https://mengxu98.github.io/scop/reference/MetabolismPlot.md)
  : Plots for metabolism pathway scoring
- [`RunscFEA()`](https://mengxu98.github.io/scop/reference/RunscFEA.md)
  : Run scFEA flux estimation for a Seurat object
- [`scFEAHeatmap()`](https://mengxu98.github.io/scop/reference/scFEAHeatmap.md)
  : Plot scFEA module flux heatmap
- [`scFEAVolcanoPlot()`](https://mengxu98.github.io/scop/reference/scFEAVolcanoPlot.md)
  : Plot scFEA flux Cohen's d volcano plots
- [`scFEABalanceBarPlot()`](https://mengxu98.github.io/scop/reference/scFEABalanceBarPlot.md)
  : Plot scFEA metabolite balance changes

## Trajectory and Cell Fate

### Trajectory and Pseudotime Inference

- [`RunMonocle2()`](https://mengxu98.github.io/scop/reference/RunMonocle2.md)
  : Run Monocle2 analysis
- [`RunMonocle3()`](https://mengxu98.github.io/scop/reference/RunMonocle3.md)
  : Run Monocle3 analysis
- [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md) :
  Run PAGA analysis
- [`PAGAPlot()`](https://mengxu98.github.io/scop/reference/PAGAPlot.md)
  : PAGA plot
- [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md)
  : Run Palantir analysis
- [`PalantirTrajectoryPlot()`](https://mengxu98.github.io/scop/reference/PalantirTrajectoryPlot.md)
  : Plot Palantir trajectories
- [`RunSlingshot()`](https://mengxu98.github.io/scop/reference/RunSlingshot.md)
  : RunSlingshot

### RNA Velocity and Vector Fields

- [`RunCell2fate()`](https://mengxu98.github.io/scop/reference/RunCell2fate.md)
  : Run Cell2fate RNA velocity analysis
- [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md)
  : Run scVelo workflow
- [`VelocityPlot()`](https://mengxu98.github.io/scop/reference/VelocityPlot.md)
  : Velocity Plot
- [`RunVECTOR()`](https://mengxu98.github.io/scop/reference/RunVECTOR.md)
  : Run VECTOR developmental direction inference
- [`VECTORPlot()`](https://mengxu98.github.io/scop/reference/VECTORPlot.md)
  : Plot VECTOR results

### Fate Mapping and Lineages

- [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md)
  : Run CellRank analysis
- [`CellRankPlot()`](https://mengxu98.github.io/scop/reference/CellRankPlot.md)
  : Plot stored CellRank results
- [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md) :
  Run WOT analysis
- [`LineagePlot()`](https://mengxu98.github.io/scop/reference/LineagePlot.md)
  : Lineage Plot
- [`BranchStreamPlot()`](https://mengxu98.github.io/scop/reference/BranchStreamPlot.md)
  : Branch Stream Plot
- [`PseudotimeProjectionPlot()`](https://mengxu98.github.io/scop/reference/PseudotimeProjectionPlot.md)
  : Pseudotime Projection Plot

### Dynamic Expression and Pathways

- [`RunDynamicFeatures()`](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md)
  : Calculates dynamic features for lineages
- [`DynamicPlot()`](https://mengxu98.github.io/scop/reference/DynamicPlot.md)
  : Plot dynamic features across pseudotime
- [`DynamicHeatmap()`](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md)
  : Heatmap plot for dynamic features along lineages
- [`RunDynamicEnrichment()`](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)
  : RunDynamicEnrichment
- [`RunCellRankTrends()`](https://mengxu98.github.io/scop/reference/RunCellRankTrends.md)
  : Fit and cluster CellRank lineage trends
- [`RunCellRankEnrichment()`](https://mengxu98.github.io/scop/reference/RunCellRankEnrichment.md)
  : Enrich CellRank trend modules

### Cell Potency and Aging

- [`PrepareCytoTRACE2()`](https://mengxu98.github.io/scop/reference/PrepareCytoTRACE2.md)
  : Prepare CytoTRACE2 model resources
- [`RunCytoTRACE()`](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md)
  : Run CytoTRACE 2
- [`CytoTRACEPlot()`](https://mengxu98.github.io/scop/reference/CytoTRACEPlot.md)
  : Plot CytoTRACE 2 Results
- [`RunFitDevo()`](https://mengxu98.github.io/scop/reference/RunFitDevo.md)
  : Run FitDevo developmental potential scoring
- [`FitDevoPlot()`](https://mengxu98.github.io/scop/reference/FitDevoPlot.md)
  : Plot FitDevo results
- [`RuntAge()`](https://mengxu98.github.io/scop/reference/RuntAge.md) :
  Run tAge transcriptomic aging-clock prediction
- [`tAgePlot()`](https://mengxu98.github.io/scop/reference/tAgePlot.md)
  : Plot tAge transcriptomic aging-clock predictions

## Cell-Cell Communication

### Communication Inference

- [`RunCCC()`](https://mengxu98.github.io/scop/reference/RunCCC.md) :
  Run common cell-cell communication analyses
- [`RunCellChat()`](https://mengxu98.github.io/scop/reference/RunCellChat.md)
  : Run CellChat analysis
- [`RunCellphoneDB()`](https://mengxu98.github.io/scop/reference/RunCellphoneDB.md)
  : Run CellphoneDB analysis
- [`RunLIANA()`](https://mengxu98.github.io/scop/reference/RunLIANA.md)
  : Run LIANA cell-cell communication analysis
- [`RunMDIC3()`](https://mengxu98.github.io/scop/reference/RunMDIC3.md)
  : Run MDIC3 cell-cell communication inference
- [`RunscOMM()`](https://mengxu98.github.io/scop/reference/RunscOMM.md)
  : Run scOMM label prediction
- [`RunNichenetr()`](https://mengxu98.github.io/scop/reference/RunNichenetr.md)
  : Run NicheNet analysis
- [`RunMultiNichenetr()`](https://mengxu98.github.io/scop/reference/RunMultiNichenetr.md)
  : Run MultiNicheNet analysis

### Spatial Communication

- [`RunSpaTalk()`](https://mengxu98.github.io/scop/reference/RunSpaTalk.md)
  : Run official SpaTalk spatial communication analysis
- [`SpaTalkPlot()`](https://mengxu98.github.io/scop/reference/SpaTalkPlot.md)
  : Plot stored SpaTalk results
- [`RunCOMMOT()`](https://mengxu98.github.io/scop/reference/RunCOMMOT.md)
  : Run official COMMOT spatial communication analysis
- [`COMMOTPlot()`](https://mengxu98.github.io/scop/reference/COMMOTPlot.md)
  : Plot stored COMMOT results
- [`RunSpatialDM()`](https://mengxu98.github.io/scop/reference/RunSpatialDM.md)
  : Run official SpatialDM spatial ligand-receptor association analysis
- [`SpatialDMPlot()`](https://mengxu98.github.io/scop/reference/SpatialDMPlot.md)
  : Plot stored SpatialDM results
- [`GetSpatialDMResult()`](https://mengxu98.github.io/scop/reference/GetSpatialDMResult.md)
  : Access stored SpatialDM results
- [`RunSpatialCellChat()`](https://mengxu98.github.io/scop/reference/RunSpatialCellChat.md)
  : Run Spatial CellChat analysis
- [`SpatialCellChatPlot()`](https://mengxu98.github.io/scop/reference/SpatialCellChatPlot.md)
  : Plot stored SpatialCellChat results

### Signaling Activity and Dynamics

- [`RunSecAct()`](https://mengxu98.github.io/scop/reference/RunSecAct.md)
  : Run SecAct secreted protein activity inference
- [`RunSecActCCC()`](https://mengxu98.github.io/scop/reference/RunSecActCCC.md)
  : Run SecAct cell-cell communication analysis
- [`RunSecActSignalingPattern()`](https://mengxu98.github.io/scop/reference/RunSecActSignalingPattern.md)
  : Run SecAct spatial signaling pattern analysis
- [`RunSecActPatternGenes()`](https://mengxu98.github.io/scop/reference/RunSecActPatternGenes.md)
  : Extract SecAct pattern-associated secreted proteins
- [`RunSecActVelocity()`](https://mengxu98.github.io/scop/reference/RunSecActVelocity.md)
  : Run SecAct spatial signaling velocity

### Communication Resources and Interoperability

- [`GetCCCObject()`](https://mengxu98.github.io/scop/reference/GetCCCObject.md)
  : Get an upstream CellChat-family object
- [`ListCCCDB()`](https://mengxu98.github.io/scop/reference/ListCCCDB.md)
  : List CCC ligand-receptor and prior-model databases
- [`PrepareCCCDB()`](https://mengxu98.github.io/scop/reference/PrepareCCCDB.md)
  : Prepare cell-cell communication databases
- [`ccc_to_adata()`](https://mengxu98.github.io/scop/reference/ccc_to_adata.md)
  : Convert CCC results to OmicVerse communication AnnData
- [`ccc_to_liana()`](https://mengxu98.github.io/scop/reference/ccc_to_liana.md)
  : Convert CCC results to a LIANA-like table

### Communication Visualization

- [`CCCHeatmap()`](https://mengxu98.github.io/scop/reference/CCCHeatmap.md)
  : CCC heatmap and dot matrix plot
- [`CCCNetworkPlot()`](https://mengxu98.github.io/scop/reference/CCCNetworkPlot.md)
  : CCC network and flow plots
- [`CCCStatPlot()`](https://mengxu98.github.io/scop/reference/CCCStatPlot.md)
  : CCC statistical distribution and summary plots

## Gene Regulatory Networks

### Network Inference

- [`RunGRN()`](https://mengxu98.github.io/scop/reference/RunGRN.md) :
  Infer gene regulatory networks with a selected backend method
- [`RunGRNBoost2()`](https://mengxu98.github.io/scop/reference/RunGRNBoost2.md)
  : Infer gene regulatory networks with GRNBoost2
- [`RunGENIE3()`](https://mengxu98.github.io/scop/reference/RunGENIE3.md)
  : Infer gene regulatory networks with GENIE3
- [`RunGNIPLR()`](https://mengxu98.github.io/scop/reference/RunGNIPLR.md)
  : Infer gene regulatory networks with GNIPLR

### Regulon Inference and Activity

- [`RunSCENIC()`](https://mengxu98.github.io/scop/reference/RunSCENIC.md)
  : Run SCENIC gene regulatory network analysis
- [`SCENICPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlot.md)
  : Plot SCENIC or SCENIC+ regulon activity
- [`RunCisTarget()`](https://mengxu98.github.io/scop/reference/RunCisTarget.md)
  : Run cisTarget motif enrichment on a GRN adjacency table
- [`RunDorothea()`](https://mengxu98.github.io/scop/reference/RunDorothea.md)
  : Run transcription factor activity inference
- [`DorotheaPlot()`](https://mengxu98.github.io/scop/reference/DorotheaPlot.md)
  : Plot transcription factor activity

### Multiomic Regulatory Networks

- [`RunSCENICPlus()`](https://mengxu98.github.io/scop/reference/RunSCENICPlus.md)
  : Read SCENIC+ results or run an approximate SCENICPlus-compatible
  eGRN analysis
- [`SCENICPlusPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlusPlot.md)
  : Plot SCENIC+ eRegulon activity

## Phenotype and Perturbation Analysis

### Phenotype Association

- [`RunScissor()`](https://mengxu98.github.io/scop/reference/RunScissor.md)
  : Run Scissor phenotype-associated cell selection
- [`ScissorPlot()`](https://mengxu98.github.io/scop/reference/ScissorPlot.md)
  : Plot Scissor results
- [`RunAugur()`](https://mengxu98.github.io/scop/reference/RunAugur.md)
  : Prioritize perturbed cell types using Augur
- [`RunFWP()`](https://mengxu98.github.io/scop/reference/RunFWP.md) :
  Run FWP feature-weight phenotype scoring

### Network Perturbation

- [`RunscTenifoldKnk()`](https://mengxu98.github.io/scop/reference/RunScTenifoldKnk.md)
  : Run scTenifoldKnk in-silico knockout analysis
- [`scTenifoldKnkPlot()`](https://mengxu98.github.io/scop/reference/scTenifoldKnkPlot.md)
  : scTenifoldKnk Plot
- [`RunscTenifoldNet()`](https://mengxu98.github.io/scop/reference/RunscTenifoldNet.md)
  : Run scTenifoldNet network comparison
- [`scTenifoldNetPlot()`](https://mengxu98.github.io/scop/reference/scTenifoldNetPlot.md)
  : scTenifoldNet Plot

## ATAC and Multiome Analysis

### Accessibility Quality and Coverage

- [`RunATACQC()`](https://mengxu98.github.io/scop/reference/RunATACQC.md)
  : Run scATAC quality control metrics
- [`CoverageTrackPlot()`](https://mengxu98.github.io/scop/reference/CoverageTrackPlot.md)
  : Coverage track plot for ATAC data

### Peak Definition and Annotation

- [`CombineATACPeaks()`](https://mengxu98.github.io/scop/reference/CombineATACPeaks.md)
  : Combine explicit peak sets into one common interval universe
- [`RecountATAC()`](https://mengxu98.github.io/scop/reference/RecountATAC.md)
  : Recount fragments over an explicit common ATAC peak universe
- [`FilterATACPeaks()`](https://mengxu98.github.io/scop/reference/FilterATACPeaks.md)
  : Filter an ATAC peak universe with explicit criteria
- [`AnnotateATACPeaks()`](https://mengxu98.github.io/scop/reference/AnnotateATACPeaks.md)
  : Annotate ATAC peaks with nearest genomic features

### Motifs and Accessibility Deviations

- [`AddATACMotifs()`](https://mengxu98.github.io/scop/reference/AddATACMotifs.md)
  : Add native motif matches and sequence statistics
- [`RunMotifEnrichment()`](https://mengxu98.github.io/scop/reference/RunMotifEnrichment.md)
  : Test motif enrichment within an explicit peak universe
- [`RunChromVAR()`](https://mengxu98.github.io/scop/reference/RunChromVAR.md)
  : Compute native chromVAR deviations from an explicit sequence
  reference

### Gene Activity and Peak-Gene Associations

- [`RunGeneActivity()`](https://mengxu98.github.io/scop/reference/RunGeneActivity.md)
  : Quantify gene activity from indexed ATAC fragments
- [`RunLinkPeaks()`](https://mengxu98.github.io/scop/reference/RunLinkPeaks.md)
  : Link ATAC peaks to paired measured expression

## Spatial Omics Analysis

### Spatial Quality Control and Normalization

- [`RunSpotQC()`](https://mengxu98.github.io/scop/reference/RunSpotQC.md)
  : Run spot-level quality control
- [`RunSpotSweeper()`](https://mengxu98.github.io/scop/reference/RunSpotSweeper.md)
  : Run SpotSweeper spatial quality control
- [`RunSpaNorm()`](https://mengxu98.github.io/scop/reference/RunSpaNorm.md)
  : Run SpaNorm spatial normalization

### Deconvolution and Ecotypes

- [`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md) :
  Run RCTD spatial deconvolution
- [`RunCell2location()`](https://mengxu98.github.io/scop/reference/RunCell2location.md)
  : Run cell2location spatial deconvolution
- [`Cell2locationPlot()`](https://mengxu98.github.io/scop/reference/Cell2locationPlot.md)
  : Plot cell2location spatial results
- [`RunCARD()`](https://mengxu98.github.io/scop/reference/RunCARD.md) :
  Run CARD spatial deconvolution
- [`RunSPOTlight()`](https://mengxu98.github.io/scop/reference/RunSPOTlight.md)
  : Run SPOTlight spatial deconvolution
- [`RunSTdeconvolve()`](https://mengxu98.github.io/scop/reference/RunSTdeconvolve.md)
  : Run STdeconvolve reference-free spatial deconvolution
- [`STdeconvolvePlot()`](https://mengxu98.github.io/scop/reference/STdeconvolvePlot.md)
  : Plot STdeconvolve topic proportions
- [`RunSpatialEcoTyper()`](https://mengxu98.github.io/scop/reference/RunSpatialEcoTyper.md)
  : Run SpatialEcoTyper spatial ecotype analysis
- [`SpatialDeconvolutionPlot()`](https://mengxu98.github.io/scop/reference/SpatialDeconvolutionPlot.md)
  : Plot stored spatial deconvolution proportions

### Spatial Domains and Cell Mapping

- [`RunBayesSpace()`](https://mengxu98.github.io/scop/reference/RunBayesSpace.md)
  : Run BayesSpace spatial clustering
- [`RunBANKSY()`](https://mengxu98.github.io/scop/reference/RunBANKSY.md)
  : Run BANKSY spatial clustering
- [`RunCytoSPACE()`](https://mengxu98.github.io/scop/reference/RunCytoSPACE.md)
  : Run CytoSPACE spatial assignment
- [`RunSmoothClust()`](https://mengxu98.github.io/scop/reference/RunSmoothClust.md)
  : Run smoothclust spatial domain clustering

### Spatial Feature Patterns

- [`RunMERINGUE()`](https://mengxu98.github.io/scop/reference/RunMERINGUE.md)
  : Run MERINGUE spatial autocorrelation analysis
- [`RunSpatialVariableFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialVariableFeatures.md)
  : Run spatial variable feature detection
- [`SpatialVariableFeaturePlot()`](https://mengxu98.github.io/scop/reference/SpatialVariableFeaturePlot.md)
  : Plot spatial variable feature results
- [`RunSpatialGradientFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialGradientFeatures.md)
  : Run spatial gradient feature screening
- [`SpatialGradientPlot()`](https://mengxu98.github.io/scop/reference/SpatialGradientPlot.md)
  : Plot spatial gradient screening results
- [`RunCSIDE()`](https://mengxu98.github.io/scop/reference/RunCSIDE.md)
  : Run C-SIDE spatial differential expression

### Spatial Neighborhoods and Interactions

- [`RunSpatialNeighborhood()`](https://mengxu98.github.io/scop/reference/RunSpatialNeighborhood.md)
  : Run spatial neighborhood statistics
- [`SpatialNeighborhoodProfile()`](https://mengxu98.github.io/scop/reference/SpatialNeighborhoodProfile.md)
  : Profile cell neighborhoods at several distances
- [`SpatialNeighborhoodPlot()`](https://mengxu98.github.io/scop/reference/SpatialNeighborhoodPlot.md)
  : Spatial neighborhood plot
- [`RunSpatialNetwork()`](https://mengxu98.github.io/scop/reference/RunSpatialNetwork.md)
  : Build a native spatial network
- [`SpatialNetworkPlot()`](https://mengxu98.github.io/scop/reference/SpatialNetworkPlot.md)
  : Plot a native spatial network
- [`RunMistyR()`](https://mengxu98.github.io/scop/reference/RunMistyR.md)
  : Run mistyR multiview spatial modeling
- [`MistyRPlot()`](https://mengxu98.github.io/scop/reference/MistyRPlot.md)
  : Plot stored MISTy results
- [`RunStatialKontextual()`](https://mengxu98.github.io/scop/reference/RunStatialKontextual.md)
  : Run Statial Kontextual spatial relationships
- [`StatialKontextualPlot()`](https://mengxu98.github.io/scop/reference/StatialKontextualPlot.md)
  : Plot stored Statial Kontextual results

### Spatial Integration and Evaluation

- [`RunSpatialIntegration()`](https://mengxu98.github.io/scop/reference/RunSpatialIntegration.md)
  : Run multi-sample spatial integration
- [`SpatialIntegrationPlot()`](https://mengxu98.github.io/scop/reference/SpatialIntegrationPlot.md)
  : Plot spatial integration results
- [`RunSpatialBenchmark()`](https://mengxu98.github.io/scop/reference/RunSpatialBenchmark.md)
  : Benchmark spatial domain clustering methods
- [`SpatialBenchmarkPlot()`](https://mengxu98.github.io/scop/reference/SpatialBenchmarkPlot.md)
  : Plot spatial domain clustering benchmarks

### Spatial Maps

- [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md)
  : Spatial spot plot
- [`SpatialCellPlot()`](https://mengxu98.github.io/scop/reference/SpatialCellPlot.md)
  : Plot spatial cell boundaries

## Genomic Variation and Trait Association

### Copy Number Variation

- [`RunCNV()`](https://mengxu98.github.io/scop/reference/RunCNV.md) :
  Run copy-number alteration inference
- [`CNVPlot()`](https://mengxu98.github.io/scop/reference/CNVPlot.md) :
  Plot SCOP CNV results

### GWAS and Trait Association

- [`RunscPagwas()`](https://mengxu98.github.io/scop/reference/RunscPagwas.md)
  : Run scPagwas
- [`PlotscPagwas()`](https://mengxu98.github.io/scop/reference/PlotscPagwas.md)
  : Plot scPagwas Scores

## Visualization and Interactive Exploration

### Feature and Statistical Plots

- [`CellStatPlot()`](https://mengxu98.github.io/scop/reference/CellStatPlot.md)
  : Statistical plot of cells
- [`FeatureCorPlot()`](https://mengxu98.github.io/scop/reference/FeatureCorPlot.md)
  : Features correlation plot
- [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)
  : Feature statistical plots

### Heatmaps and Themes

- [`CellCorHeatmap()`](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md)
  : Cell Correlation Heatmap
- [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md)
  : Feature Heatmap
- [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)
  : Group heatmap
- [`theme_scop()`](https://mengxu98.github.io/scop/reference/theme_scop.md)
  : Default ggplot2 theme for scop plots
- [`theme_spatial()`](https://mengxu98.github.io/scop/reference/theme_spatial.md)
  : ggplot2 theme for spatial plots

### Interactive Exploration

- [`PrepareSCExplorer()`](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md)
  : Prepare Seurat objects for the SCExplorer
- [`RunSCExplorer()`](https://mengxu98.github.io/scop/reference/RunSCExplorer.md)
  : Run SCExplorer

## Data I/O and Object Management

### Data Import and Validation

- [`ReadATACData()`](https://mengxu98.github.io/scop/reference/ReadATACData.md)
  : Read raw ATAC counts into a ChromatinAssay

- [`ReadMultiomeData()`](https://mengxu98.github.io/scop/reference/ReadMultiomeData.md)
  : Read paired 10x RNA and ATAC counts

- [`ReadSpatialData()`](https://mengxu98.github.io/scop/reference/ReadSpatialData.md)
  : Load and validate Visium, Visium HD, or Xenium data

- [`CheckATACData()`](https://mengxu98.github.io/scop/reference/CheckATACData.md)
  : Check an ATAC object without changing it

- [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  :

  Check and report the type of data in `Seurat` object

- [`CheckDataList()`](https://mengxu98.github.io/scop/reference/CheckDataList.md)
  :

  Check and preprocess a list of `Seurat` objects

- [`CheckDataMerge()`](https://mengxu98.github.io/scop/reference/CheckDataMerge.md)
  : Check and preprocess a merged seurat object

- [`ExportATACResults()`](https://mengxu98.github.io/scop/reference/ExportATACResults.md)
  : Export an ATAC object and explicit result/resource tables

### Object and File Conversion

- [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md)
  : Convert an anndata object to a seurat object

- [`srt_to_adata()`](https://mengxu98.github.io/scop/reference/srt_to_adata.md)
  : Convert a Seurat object to an AnnData object

- [`h5ad_to_srt()`](https://mengxu98.github.io/scop/reference/h5ad_to_srt.md)
  :

  Read an `.h5ad` file and convert to a `Seurat`

- [`srt_to_h5ad()`](https://mengxu98.github.io/scop/reference/srt_to_h5ad.md)
  :

  Convert a Seurat object to an `.h5ad` file

- [`loom_to_adata()`](https://mengxu98.github.io/scop/reference/loom_to_adata.md)
  :

  Read a `.loom` file as an AnnData object

- [`loom_to_srt()`](https://mengxu98.github.io/scop/reference/loom_to_srt.md)
  :

  Read a `.loom` file and convert to a `Seurat`

### Spatial Framework Interoperability

- [`srt_to_spe()`](https://mengxu98.github.io/scop/reference/srt_to_spe.md)
  : Convert Seurat to SpatialExperiment
- [`spe_to_srt()`](https://mengxu98.github.io/scop/reference/spe_to_srt.md)
  : Convert SpatialExperiment to Seurat
- [`srt_to_giotto()`](https://mengxu98.github.io/scop/reference/srt_to_giotto.md)
  : Convert Seurat to a native Giotto object
- [`giotto_to_srt()`](https://mengxu98.github.io/scop/reference/giotto_to_srt.md)
  : Convert Giotto to Seurat

### Spatial Coordinates and Results

- [`SpatialCoordinates()`](https://mengxu98.github.io/scop/reference/SpatialCoordinates.md)
  : Read spatial coordinates with an explicit coordinate contract
- [`SetSpatialImageAxes()`](https://mengxu98.github.io/scop/reference/SetSpatialImageAxes.md)
  : Persist the coordinate axis convention of a VisiumV2 image
- [`GetSpatialGraph()`](https://mengxu98.github.io/scop/reference/GetSpatialGraph.md)
  : Read or convert a stored spatial graph
- [`GetSpatialResult()`](https://mengxu98.github.io/scop/reference/GetSpatialResult.md)
  : Access spatial analysis results

### Object Access and Organization

- [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
  : Find the default reduction name in a Seurat object

- [`FetchDataZero()`](https://mengxu98.github.io/scop/reference/FetchDataZero.md)
  : FetchData but with zeroes for unavailable genes

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

### Feature Annotation and Identifiers

- [`AddFeaturesData()`](https://mengxu98.github.io/scop/reference/AddFeaturesData.md)
  : Add features data
- [`GetFeaturesData()`](https://mengxu98.github.io/scop/reference/GetFeaturesData.md)
  : Get features data
- [`AnnotateFeatures()`](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md)
  : Annotate Features
- [`GetSimilarFeatures()`](https://mengxu98.github.io/scop/reference/GetSimilarFeatures.md)
  : Find features with expression patterns similar to provided features
- [`RenameFeatures()`](https://mengxu98.github.io/scop/reference/RenameFeatures.md)
  : Rename features for the Seurat object
- [`GeneConvert()`](https://mengxu98.github.io/scop/reference/GeneConvert.md)
  : Gene ID conversion function using biomart
- [`ConvertHomologs()`](https://mengxu98.github.io/scop/reference/ConvertHomologs.md)
  : Convert homologous gene symbols in expression objects

### Interactive Data Files

- [`CreateDataFile()`](https://mengxu98.github.io/scop/reference/CreateDataFile.md)
  : Create HDF5 data file from Seurat object
- [`CreateMetaFile()`](https://mengxu98.github.io/scop/reference/CreateMetaFile.md)
  : Create Meta File in HDF5 format from Seurat object
- [`FetchH5()`](https://mengxu98.github.io/scop/reference/FetchH5.md) :
  Fetch data from the hdf5 file and returns a Seurat object

## Resources and Environment

### Resource Preparation

- [`ListDB()`](https://mengxu98.github.io/scop/reference/ListDB.md) :
  List cached gene annotation databases
- [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md)
  : Prepare databases and reference resources
- [`PrepareCustomDB()`](https://mengxu98.github.io/scop/reference/PrepareCustomDB.md)
  : Prepare custom gene annotation databases

### Pathway and Gene-Set Resources

- [`PrepareGO()`](https://mengxu98.github.io/scop/reference/PrepareGO.md)
  : Prepare GO databases
- [`PrepareKEGG()`](https://mengxu98.github.io/scop/reference/PrepareKEGG.md)
  : Prepare KEGG databases
- [`PrepareWikiPathway()`](https://mengxu98.github.io/scop/reference/PrepareWikiPathway.md)
  : Prepare WikiPathway databases
- [`PrepareReactome()`](https://mengxu98.github.io/scop/reference/PrepareReactome.md)
  : Prepare Reactome databases
- [`PrepareMSigDB()`](https://mengxu98.github.io/scop/reference/PrepareMSigDB.md)
  : Prepare MSigDB databases
- [`PrepareCORUM()`](https://mengxu98.github.io/scop/reference/PrepareCORUM.md)
  : Prepare CORUM databases

### Annotation and Ontology Resources

- [`PrepareMP()`](https://mengxu98.github.io/scop/reference/PrepareMP.md)
  : Prepare MP databases
- [`PrepareDO()`](https://mengxu98.github.io/scop/reference/PrepareDO.md)
  : Prepare DO databases
- [`PrepareHPO()`](https://mengxu98.github.io/scop/reference/PrepareHPO.md)
  : Prepare HPO databases
- [`PreparePFAM()`](https://mengxu98.github.io/scop/reference/PreparePFAM.md)
  : Prepare PFAM databases
- [`PrepareChromosome()`](https://mengxu98.github.io/scop/reference/PrepareChromosome.md)
  : Prepare Chromosome databases
- [`PrepareGeneType()`](https://mengxu98.github.io/scop/reference/PrepareGeneType.md)
  : Prepare GeneType databases
- [`PrepareEnzyme()`](https://mengxu98.github.io/scop/reference/PrepareEnzyme.md)
  : Prepare Enzyme databases

### Regulatory Resources

- [`PrepareTF()`](https://mengxu98.github.io/scop/reference/PrepareTF.md)
  : Prepare TF databases
- [`PrepareTFLink()`](https://mengxu98.github.io/scop/reference/PrepareTFLink.md)
  : Prepare TFLink databases
- [`PrepareHTFtarget()`](https://mengxu98.github.io/scop/reference/PrepareHTFtarget.md)
  : Prepare hTFtarget databases
- [`PrepareTRRUST()`](https://mengxu98.github.io/scop/reference/PrepareTRRUST.md)
  : Prepare TRRUST databases
- [`PrepareJASPAR()`](https://mengxu98.github.io/scop/reference/PrepareJASPAR.md)
  : Prepare JASPAR databases
- [`PrepareENCODE()`](https://mengxu98.github.io/scop/reference/PrepareENCODE.md)
  : Prepare ENCODE databases

### Cell Surface and Signaling Resources

- [`PrepareCSPA()`](https://mengxu98.github.io/scop/reference/PrepareCSPA.md)
  : Prepare CSPA databases
- [`PrepareSurfaceome()`](https://mengxu98.github.io/scop/reference/PrepareSurfaceome.md)
  : Prepare Surfaceome databases
- [`PrepareSPRomeDB()`](https://mengxu98.github.io/scop/reference/PrepareSPRomeDB.md)
  : Prepare SPRomeDB databases
- [`PrepareVerSeDa()`](https://mengxu98.github.io/scop/reference/PrepareVerSeDa.md)
  : Prepare VerSeDa databases
- [`PrepareCellTalk()`](https://mengxu98.github.io/scop/reference/PrepareCellTalk.md)
  : Prepare CellTalk databases
- [`PrepareCellChat()`](https://mengxu98.github.io/scop/reference/PrepareCellChat.md)
  : Prepare CellChat databases

### Example Datasets

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

### Python Environment Management

- [`check_python()`](https://mengxu98.github.io/scop/reference/check_python.md)
  : Check and install python packages
- [`env_requirements()`](https://mengxu98.github.io/scop/reference/env_requirements.md)
  : Python environment requirements
- [`env_info()`](https://mengxu98.github.io/scop/reference/env_info.md)
  : Print environment information
- [`ListEnv()`](https://mengxu98.github.io/scop/reference/ListEnv.md) :
  List conda-compatible Python environments
- [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
  : Prepare the python environment
- [`remove_python()`](https://mengxu98.github.io/scop/reference/remove_python.md)
  : Remove Python packages from a conda-compatible Python environment
- [`RemoveEnv()`](https://mengxu98.github.io/scop/reference/RemoveEnv.md)
  : Remove a conda-compatible Python environment
