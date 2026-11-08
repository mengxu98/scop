# ambiguous annotation and duplicate enrichment fail clearly

    Code
      cellstate_fixture_plot(f, class.by = "annotation")
    Condition
      Error in `FUN()`:
      ! Mixed annotations in cluster A; supply an explicit map.

---

    Code
      cellstate_fixture_plot(f)
    Condition
      Error in `cellstate_gsea()`:
      ! Duplicate GSEA cluster/pathway rows; select one result per combination.
