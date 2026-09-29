# IREA experimental reproduction in SCOP

The Immune Dictionary provides mouse lymph-node single-cell profiles after 86
cytokine perturbations, plus mouse and human-orthologue signature spreadsheets.
The reference files come from <https://www.immune-dictionary.org/download/>.
The original methods are described in Cui et al., *Nature* 625:377–384 (2024),
<https://doi.org/10.1038/s41586-023-06816-9>.

`RunIREA()` currently implements the published analysis outline, but **is not
numerically equivalent to the web portal**. In particular, the portal's
reference-cell selection, hypergeometric background, and matrix filtering are
not fully specified by the paper or downloads.
Every result carries `validation = "method_reconstructed; portal numeric
concordance not established"`. Do not describe these results as a verified
official IREA run.

Run `Rscript inst/irea/irea_download.R D:/scop-dev/datasets/IREA` from the
checkout to cache the public resources and record source URLs, file sizes and
MD5 checksums. The cache is intentionally outside the SCOP package. Then run
`Rscript inst/irea/irea_demo.R D:/scop-dev/datasets/IREA output-directory` to
produce NK gene-list and matrix example tables and PDF figures.

The reference loader supports the 14 cell types that have downloadable Seurat
objects. `RunIREA()` accepts a gene list, a named contrast vector or a Seurat
object with case and control groups. The human mode uses the published
`Gene_Human` columns to map human symbols back to mouse reference genes; genes
without an unambiguous mapping are omitted and this partial map is returned.
Neither the mapping nor its numerical results have been matched to the portal.
All 14 downloaded objects were opened and produced 86 cytokine-response rows
with a matching mouse gene; each human-orthologue table also completed one
mapped-gene smoke test. These are runtime checks, not numerical replication.

The comparison uses the nine mouse example genes shown on the portal
and both columns of its `gene_test.xlsx`. The cached portal and local tables,
figures, and `numeric_comparison.csv` are in
`D:/scop-dev/datasets/IREA`. The repaired comparator matches all 360 terms,
including Greek-letter aliases, and checks effect, raw P value and FDR
separately. After the PBS correction, 81 rows agree on all three statistics;
most are zero-overlap hypergeometric entries. Of these, four nonzero rows are
the NK gene-list score polarization example. The other modes still do not
pass numerical concordance. As a concrete example,
the portal gives the IL-15 NK gene-list score as 1.5865844168, while this
implementation gives 1.8905808575; the portal FDR is 3.25e-5 versus local
4.53e-15. For the anti-PD-1 matrix example, portal IL-18 effect is
0.1231831905 versus local 0.0898199336. These differences prevent a claim of
complete reproduction. Check current values in the CSV rather than treating
the numbers in this note as a permanent baseline.

## September 2026 repair and reproducible validation

Polarization score/projection previously compared the target state against
cells labelled `None`. They now use PBS-treated cells as the control. In the
cached nine-gene NK score example, this recovers all four portal effects,
P values and FDR values to the exported precision. This evidence covers that
mouse NK gene-list example only: it does not establish equivalence for human
mapping, other cell types, hypergeometric analysis or matrix projection.

The Wilcoxon helper requests both one-sided asymptotic tail probabilities and
returns twice the smaller tail, capped at one, with continuity/tie corrections
and explicit full-precision ranking. This avoids cancellation in R 4.6.1's
two-sided `1 - pnorm(z)` path, which otherwise turns the NK-a fixture's P value
of approximately 4.88e-208 into zero. The strict portal tolerances are retained.

The matrix portal tables also pool P values across the two example contrasts
when applying BH correction. Matrix files/data frames now default to joint
adjustment across all supplied columns before returning the requested
`contrast`; set `fdr_scope="selected_contrast"` to request per-contrast BH.
The actual family and adjusted contrast names are saved in `parameters`.
This fixes the adjustment family but does not resolve the projection/P-value
differences, so matrix FDR values remain unvalidated against the portal.

Seurat input now aligns grouping metadata to the selected layer's cell names,
ignores cells with missing group labels, and requires both groups to exist in
that layer. An exact layer name is required: split layers must be joined before
comparing across batches. The function does not silently use the first layer.
Different scoring methods or species cannot share a comparison plot. Radar
scores are recomputed using the plotting `fdr_cutoff` without changing the
saved statistical table. Plots identify the implementation as experimental.

The review outputs are in `D:/scop-dev/datasets/IREA/repair_20260929`.
The original cached results are retained for before/after comparison. Run:

```sh
python inst/irea/irea_compare.py D:/scop-dev/datasets/IREA --local-dir D:/scop-dev/datasets/IREA/repair_20260929
python inst/irea/irea_compare.py D:/scop-dev/datasets/IREA --local-dir D:/scop-dev/datasets/IREA/repair_20260929 --analysis cell_polarization --method score --output D:/scop-dev/datasets/IREA/repair_20260929/polarization_score_concordance.csv --require-concordance
```

The second command must pass for this fixed NK fixture. Applying
`--require-concordance` to all modes must currently fail. A match of an effect
alone, or a small absolute difference between extremely small P values, does
not count as full concordance. P values/FDR use relative tolerance of 1e-8
plus rounding tolerance derived from the last printed digit of older cached
tables, rather than a fixed absolute tolerance. Effects use absolute and
relative tolerance 1e-8. New portal captures preserve full-precision HTML
sort values when available.

Two new submissions of the public nine-gene response example on 2026-09-29
returned identical results and agreed with the earlier cached effects. The
remaining response discrepancy is reproducible, rather than merely a stale
cached website result. No private user data were submitted.

The downloaded matrix example covers the anti-PD-1 and anti-CTLA-4 contrasts
for NK cells. The original tumour single-cell dataset has not yet been
reprocessed for macrophage and CD8 T-cell polarization; those biological
comparisons remain unverified.

The response score represents similarity to a perturbation signature, not a
measured cytokine concentration or proof of causal action. Human mode still
uses the original mouse perturbation data as the experimental reference.
