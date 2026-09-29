# IREA experimental reproduction in SCOP

The Immune Dictionary provides mouse lymph-node single-cell profiles after 86
cytokine perturbations, plus mouse and human-orthologue signature spreadsheets.
The reference files come from <https://www.immune-dictionary.org/download/>.
The original methods are described in Cui et al., *Nature* 625:377–384 (2024),
<https://doi.org/10.1038/s41586-023-06816-9>.

`RunIREA()` currently implements the published analysis outline, but **is not
numerically equivalent to the web portal**. In particular, the portal's
reference-cell subsampling, hypergeometric background, matrix filtering, and
polarization comparator are not fully specified by the paper or downloads.
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

The current comparison used the nine mouse example genes shown on the portal
and both columns of its `gene_test.xlsx`. The cached portal and local tables,
figures, and `numeric_comparison.csv` are in
`D:/scop-dev/datasets/IREA`. Of 360 terms, 304 could be matched by a
display-name mapping. Only 74 effects were numerically identical at 1e-8;
most were zero-overlap hypergeometric entries, and the nonzero matches are
overlap fractions with discrepant P values. As a concrete example,
the portal gives the IL-15 NK gene-list score as 1.5865844168, while this
implementation gives 1.8905808575; the portal FDR is 3.25e-5 versus local
4.53e-15. For the anti-PD-1 matrix example, portal IL-18 effect is
0.1231831905 versus local 0.0898199336. These differences prevent a claim of
complete reproduction. Check current values in the CSV rather than treating
the numbers in this note as a permanent baseline.

The downloaded matrix example covers the anti-PD-1 and anti-CTLA-4 contrasts
for NK cells. The original tumour single-cell dataset has not yet been
reprocessed for macrophage and CD8 T-cell polarization; those biological
comparisons remain unverified.

The response score represents similarity to a perturbation signature, not a
measured cytokine concentration or proof of causal action. Human mode still
uses the original mouse perturbation data as the experimental reference.
