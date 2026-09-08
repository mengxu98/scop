# Spatial platform acceptance data

Run `Rscript scripts/validate-spatial-platforms.R DATA_ROOT OUTPUT_DIR` from a
compiled SCOP checkout. It runs one isolated process per case and saves
`acceptance.csv`, result RDS files, PNGs, logs and session information. Cases may
also be named after the two directory arguments. A non-passed case returns a
nonzero exit status. Peak memory is process-tree RSS sampled every 250 ms,
including package loading; it is not a continuously measured hardware peak.

No input files are redistributed in the package. Downloads are public fixtures;
respect their original attribution and licenses. The validation runner does not
download files or install missing optional backends.

| Fixture | DATA_ROOT location | Source and interpretation |
|---|---|---|
| Visium pancreas | bundled `visium_human_pancreas_sub` | Existing SCOP dataset, GSE254829; subsets retain real expression, labels and coordinates |
| Visium HD tiny mouse brain | `hd/` vendor outputs | 10x Space Ranger 4.0.1 developer fixture, downsampled data; format/computation testing only |
| Xenium tiny mouse ileum | `xenium/` vendor outputs | 10x Xenium 3.0.0 tiny fixture, true vendor cell boundaries; not a cohort |
| Xenium human pancreas | `xenium_human_pancreas_sub.rds` | `mengxu98/datasets`, derived from TENxXeniumData; existing 3000-cell centroid subset without polygon boundaries |
| Damond IMC sample design | `diabetesData.rda` | SydneyBioX/spicyR public cell metadata; real patient/image/condition identity, no synthetic expression matrix |

## Downloads

HD documentation and license (CC BY 4.0):
<https://www.10xgenomics.com/support/software/space-ranger/latest/resources/visium-hd-example-data>

HD archive (311404970 bytes; published MD5 `7cab710801d3776125a7aa5a792cd7e4`):
<https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Tiny_3prime_Dataset/Visium_HD_Tiny_3prime_Dataset_outs.zip>

Xenium archive, also used in the Bioconductor XeniumIO import vignette:
<https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip>

Pancreas subset:
<https://github.com/mengxu98/datasets/tree/main/Xenium>

IMC data pinned to spicyR commit `64675c3eadc689267f6116017396a02717386ff3`:
<https://raw.githubusercontent.com/SydneyBioX/spicyR/64675c3eadc689267f6116017396a02717386ff3/data/diabetesData.rda>

The IMC validation uses the actual `case`, `imageID`, `stage` and `cellType`
columns. The source contains 12 patients; the selected Non-diabetic/Onset
contrast contains four patients per group. Long-duration patients are not part
of that contrast. Test results demonstrate implementation and replication
handling, not a new biological claim or an adjusted disease model.

## Acceptance scope

The runner verifies loader output identity and provenance, installed BANKSY,
smoothclust, SpaNorm and local SpotSweeper execution through the standard
workflow, native SVG, sketch/full projection, segmentation QC, subject-level
comparisons, save/reload equality and plotting. Artifact detection is explicitly
off in the SpotSweeper case. Sparse SPARK-X adapter behavior is unit-tested;
SPARK/nnSVG, arbitrary on-disk matrix backends, full-size production slides,
physical slice registration and 3D reconstruction require separate integration
validation. Do not present unavailable backends or developer fixtures as
biological benchmark evidence.
