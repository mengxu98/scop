# Spatial Data and Backend Staging Notes

This note keeps spatial example data and optional spatial method backends organized while the spatial module is expanded in small PRs.

## Data Flow

SCOP should use two data paths:

| Path | Purpose | Examples | Rule |
| --- | --- | --- | --- |
| Package data in `data/` | Small, frequently used objects that should work offline in examples and tests | `pancreas_sub`, `panc8_sub`, `pbmcmultiome_sub`, `visium_human_pancreas_sub` | Document in `R/data.R`; keep size and license acceptable for package installation. |
| External assets in `mengxu98/datasets` | Larger, derived, optional, or license-sensitive resources loaded on demand | `CytoTRACE2`, `scFEA`, `CIBERSORT`, `tAge`, `Xenium` | Store a README, manifest, checksums, size metadata, and a reproducible preparation script when possible. Cache locally under `tools::R_user_dir("scop", "data")`. |

The Xenium pancreas example should stay in `mengxu98/datasets` first. SCOP already ships a Visium pancreas object as package data, while Xenium adds an imaging-based example with separate provenance and usage constraints. Keeping it external avoids growing the base package and matches the existing on-demand pattern used by `CytoTRACE2`, `scFEA`, `CIBERSORT`, and `tAge`.

Promote an external asset to package data only when all of these are true:

- It is used by many offline examples or tests.
- The final object is small enough for normal package installation.
- The source license and redistribution terms are clear.
- The generation code is reproducible and documented in `R/data.R`.

## External Dataset Contract

Each collection in `mengxu98/datasets` should provide:

- `README.md` with source, intended SCOP use, regeneration steps, and license or usage notes.
- `manifest.tsv` with stable dataset IDs, file names, object type, source, size, checksum, and core dimensions where applicable.
- A script under `scripts/` when the asset is derived from public data.
- No empty manifest entries presented as completed datasets.

SCOP-side loaders should:

- Accept a local directory or remote base URL for reproducible tests and mirrors.
- Download only when the file is missing or fails validation.
- Validate size and checksum before reading.
- Cache under `tools::R_user_dir("scop", "data")`.
- Keep data loading separate from core method execution unless the method genuinely needs that resource.

## Spatial Backend Tiers

Spatial backends should be exposed according to implementation readiness, not according to a method wishlist.

| Tier | Meaning | Examples |
| --- | --- | --- |
| Core package data and R-native helpers | Runs with SCOP's normal package surface or small optional Bioconductor dependencies | spatial plotting helpers, coordinate parsing, `SingleCellExperiment` / `SpatialExperiment` conversion, lightweight deconvolution helpers |
| Lazy R/Bioconductor backend | Optional dependency, checked at runtime, with stable `srt@tools[[method]]` output | `CARD`, `STdeconvolve`, `BANKSY`, `BayesSpace`, `spicyR`, `Statial`, `mistyR` |
| Lazy GitHub R backend | Optional GitHub package, isolated behind a wrapper and clear dependency message | `SpatialQM` |
| Heavy bridge backend | Requires Python, graph learning, GPU, or broad multi-omics infrastructure | `SPIRAL`, `SpaMTP`, `Squidpy`, `COMMOT`, `PASTE`, `Tangram`, `cell2location` |

Do not add a method to a public `method =` argument, pkgdown reference group, or runnable example until it has:

- Runtime dependency checks with a useful install message.
- A minimal example that can be skipped safely when the optional backend is absent.
- Stable storage under `srt@tools[[method]]$parameters` and, where applicable, `srt@tools[[method]]$summary`.
- A documented output schema.
- Focused tests that use a small object or mocked backend calls.
- Plotting that follows the shared spatial visual language.

## Current Coverage

The spatial module already has or is staging coverage for the common first-pass workflows:

- QC and normalization: `SpotSweeper`, `SpaNorm`, and `SpatialQM`.
- Spatial signals: shared spatial variable feature wrappers and plotting helpers.
- Deconvolution: `RCTD`, `SPOTlight`, `CARD`, `STdeconvolve`, and `SpatialDWLS`.
- Domains: `BayesSpace`, `BANKSY`, and `SmoothClust`.
- Neighborhood and context analysis: `spicyR`, `Statial`, and `mistyR`.
- Data exchange: `SingleCellExperiment` / `SpatialExperiment` conversion and external dataset loading.

This means new wrappers should fill a clear gap rather than duplicate the same user question with a different backend.

## Next Candidate Queue

The next R-side candidates should start as smoke-test branches before they are exposed as public SCOP functions:

| Candidate | Gap filled | Gate before wrapper |
| --- | --- | --- |
| `SpaceTrooper` | Additional spatial QC workflow | Confirm stable per-cell or per-spot QC outputs and small example runtime. |
| `CRAWDAD` | Directional adjacency and neighborhood distributions | Confirm it can run from simple coordinate and label tables. |
| `PRECAST` | Multi-sample embedding, alignment, and domain analysis | Define a multi-sample output schema compatible with existing spatial integration summaries. |
| `BASS` | Shared domains across samples | Check installation/runtime cost and whether domain outputs can be summarized consistently. |
| `SpatialMNN` | Shared niche comparison across slides | Wait until SCOP's integration result schema settles. |
| `SpaTalk` / `NicheNet` | Spatial communication hypotheses | Wait for a clearer shared communication result schema. |
| `DenoIST` / `cellAdmix` | Imaging-based denoising or admixture correction | Wait for stronger Xenium/CosMx example coverage and clear correction provenance. |

## Evaluated But Deferred

`SPIRAL` is valuable for spatial representation learning, batch correction, and coordinate alignment, but it is a Python-oriented backend with additional graph-learning setup. It should wait for a clear Python dependency strategy.

`SpaMTP` is an R package, but it is a broad spatial multi-omics platform with heavy dependencies and many possible entry points. It should wait until SCOP has a narrow use case and a small wrapper target.

`HoodscanR` should remain a roadmap item until SCOP has a concrete, tested wrapper. It should not appear as a user-selectable `method =` value before that wrapper exists.

Python-side candidates such as `Squidpy`, `COMMOT`, `PASTE`, `Tangram`, and `cell2location` should remain planning items until the optional Python environment story is explicit.

## Documentation Order

Method pages come before the main tutorial:

1. Stabilize each wrapper, output schema, summary fields, examples, and dependency notes.
2. Keep planned backends in roadmap notes, not in runnable API choices.
3. Write the main spatial tutorial only after the individual method pages are stable.
