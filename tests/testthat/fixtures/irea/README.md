# NK polarization score regression fixture

`nk_polarization_score.rds` stores per-cell sums of normalized expression for
the public IREA nine-gene example, sample/polarization labels, and the four
portal output rows. It contains no private user data and no full expression
matrix. It tests numerical statistics, not just output shape or direction.

Source: Cui et al., *Dictionary of immune responses to cytokines at single-cell
resolution*, Nature 625, 377-384 (2024), DOI: 10.1038/s41586-023-06816-9.
Reference and public example: https://www.immune-dictionary.org/ .

Reference file: `downloadableData/ligands-seurat-NK_cell.RDS`.
Reference MD5: `f4eabfa734a61ad839879fd84aa7a8e9`.
Portal table captured on 2026-09-28; fixture generated on 2026-09-29.
Genes: Ifitm3, Isg15, Ifit3, Bst2, Slfn5, Isg20, Phf11b, Zbp1, Rtp4.

The score is `Matrix::colSums()` of those genes in the reference RNA data
layer, in reference cell order. Labels come from `sample` and `polarization`
metadata. Expected effects, P values and BH FDR are copied from the portal
gene-list/score/cell-polarization output. PBS is the comparator. Comparing
against `None` instead must fail the fixture. No fitting or rescaling of
expected statistics is performed.

This establishes a regression baseline only for the mouse NK nine-gene score
example. It does not validate other reference cell types or analysis modes.
