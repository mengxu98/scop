# Call MACS2 peaks for explicitly selected capture groups

Split selected cells with Signac, run the installed MACS2 executable,
and retain separate peak sets, command arguments and logs. Peak
inference is performed by MACS2; no caller is installed or silently
substituted.

## Usage

``` r
CallATACPeaks(
  object,
  sample.by,
  group.by = NULL,
  cells = NULL,
  assay = NULL,
  macs2.path = "macs2",
  outdir,
  effective.genome.size,
  qvalue = 0.05,
  extsize = 200L,
  shift = -extsize/2,
  verbose = TRUE
)
```

## Arguments

- object:

  A validated Seurat object containing a ChromatinAssay and local
  indexed fragments with explicit cell mappings.

- sample.by:

  Metadata column identifying capture/library IDs.

- group.by:

  Optional metadata column subdividing each capture, for example a cell
  type. Calls are capture-by-group, never pooled across captures.

- cells:

  Explicit selected assay cell IDs, or NULL for all assay cells.

- assay:

  Source ChromatinAssay; NULL selects the active assay.

- macs2.path:

  Path or PATH command name of an installed MACS2 2.x executable in the
  current runtime. On Windows use a prepared Linux container with
  explicit mounted paths. This function does not manage containers or
  translate host paths.

- outdir:

  New output directory. Existing paths are rejected. Split BEDs, MACS2
  outputs and logs remain here, including when a call fails.

- effective.genome.size:

  Explicit effective genome size passed to MACS2.

- qvalue:

  MACS2 q-value cutoff, strictly between zero and one.

- extsize, shift:

  MACS2 fragment extension size and shift. The supported path uses BED
  input and `--nomodel`; MACS3 FRAG arguments are not used.

- verbose:

  Whether Signac reports fragment splitting progress.

## Value

A list with named `peaks` (`GRanges` per call), `groups` with columns
`call_id`, `sample` and optional `group`, and `receipt` containing
actual cell membership, reference, fragment resources, dependency
versions, command arguments, logs, exit status and output files. The
receipt identifies the original `sample.by` and `group.by` metadata
fields; their names cannot collide with the fixed output columns. Counts
represent the supplied fragment records, not fifth-column read-support
weights. Peak significance is a caller result, not a
biological-replicate condition comparison.

## Details

Missing input, nonzero caller exit, missing output and malformed output
are errors. A successful zero-byte narrowPeak file is a valid empty peak
set; this is distinct from a failed call. No matrix or existing result
is changed. `receipt.rds` and `groups.rds` are written before
splitting/calling, and `commands.rds` plus the final `receipt.rds` are
persisted once when the function exits (including on error). These
retain selected cells, capture/group membership and resources for
diagnosing failed calls. Combine the returned sets with
[`CombineATACPeaks()`](https://mengxu98.github.io/scop/reference/CombineATACPeaks.md)
and recount the resulting ranges with
[`RecountATAC()`](https://mengxu98.github.io/scop/reference/RecountATAC.md)
before comparing previously different universes. Signac's splitting API
is reused. The process is invoked directly so its exit status and both
output streams can be retained, rather than relying on a caller wrapper
that ignores the process exit status.

## See also

[`RunATACQC()`](https://mengxu98.github.io/scop/reference/RunATACQC.md),
[`CombineATACPeaks()`](https://mengxu98.github.io/scop/reference/CombineATACPeaks.md),
[`RecountATAC()`](https://mengxu98.github.io/scop/reference/RecountATAC.md)

## Examples

``` r
if (FALSE) { # \dontrun{
called <- CallATACPeaks(atac,
  sample.by = "capture", group.by = "celltype",
  macs2.path = "/usr/local/bin/macs2", outdir = "atac_peak_calls",
  effective.genome.size = 2.7e9
)
common <- CombineATACPeaks(called$peaks, genome = "hg38")
} # }
```
