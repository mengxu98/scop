# Spatial platform acceptance data

Run `Rscript scripts/validate-spatial-platforms.R DATA_ROOT OUTPUT_DIR` from a
compiled SCOP checkout. Seven cases cover Visium baseline, BANKSY, SmoothClust,
SpaNorm, local SpotSweeper, HD import and Xenium import/polygon rendering.
Each runs in an isolated process, recording RDS, PNG, logs, session information,
elapsed time and process-tree RSS sampled every 250 ms.

The runner does not install optional backends or download data. Visium cases
use the bundled visium_human_pancreas_sub (GSE254829).

Extract the HD archive into DATA_ROOT/hd:
https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Tiny_3prime_Dataset/Visium_HD_Tiny_3prime_Dataset_outs.zip

Published MD5: 7cab710801d3776125a7aa5a792cd7e4; size: 311404970 bytes.
Source and CC BY 4.0 attribution:
https://www.10xgenomics.com/support/software/space-ranger/latest/resources/visium-hd-example-data

Extract the Xenium archive into DATA_ROOT/xenium:
https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip

Developer fixtures validate file formats and computation, not biological
performance on full production slides. HD verification checks both resolutions
and matching assays/images. Xenium uses real vendor polygons in the existing
renderer. No new segmentation QC, sketch pipeline or patient-level comparison
is included in this change.

SpotSweeper artifact detection is explicitly off in the local-QC case.
SPARK-X sparse adapter behavior is unit-tested; SPARK/nnSVG live compatibility
and arbitrary on-disk backends require separate verification. Input data are not
redistributed in SCOP. Preserve original attribution and licenses.
