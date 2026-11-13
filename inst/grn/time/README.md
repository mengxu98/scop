# Time-method implementations

`run.R` executes one validated expression and pseudotime request in an independent R process. The input matrix contains genes in rows and cells in columns. Pseudotime rows use cell IDs. LEAP accepts multiple trajectories. SCODE and SINCERITIES require one complete selected trajectory.

LEAP calls `LEAP::MAC_counter()` and orders each trajectory by time and cell ID. Its scores and signs are lagged correlations. Branch predictions retain the strongest absolute correlation for each ordered gene pair.

SCODE uses `SCODE.R` from [hmatsu1226/SCODE](https://github.com/hmatsu1226/SCODE), as supplied in the multiCSN manuscript implementation. Copyright 2016 Hirotaka Matsumoto, MIT license. Source comments were removed and `MASS::ginv()` is qualified. The numerical routine is otherwise unchanged. Each repeated run has a recorded seed. The averaged A matrix has target genes in rows and regulators in columns.

SINCERITIES uses the manuscript's parallel ridge-regression implementation derived from [CABSEL/SINCERITIES](https://github.com/CABSEL/SINCERITIES). Copyright 2016 Nan Papili Gao and Rudiyanto Gunawan, BSD 3-Clause license. The manuscript modifications are Copyright 2026 Xu Meng, MIT license. Original source notices are distributed alongside these files. The regression matrix has predictor genes in rows and response genes in columns. Edge extraction follows that orientation. This corrects the manuscript export helper's reversed endpoint labels. The complete fit is retained in `native$fit`.

SINCERITIES `SIGN = 0` leaves signs unknown. `SIGN = 1` uses the sign of partial Spearman correlation. `max_time_points = NA` retains measured times. A positive cap bins cells by their time rank after deterministic ordering. The actual times or bins and normalization divisor are stored with the result. Default edge weights match the manuscript's absolute normalized score. Signed coefficients remain available in `signed_weight`.

Algorithm errors propagate as failed runs. Valid zero scores are retained. The result does not depend on the manuscript repository at runtime.
