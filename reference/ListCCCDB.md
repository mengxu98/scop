# List CCC ligand-receptor and prior-model databases

Enumerates the ligand-receptor databases and prior models exposed by the
installed cell-cell communication methods: LIANA resources, CellChat
databases, NicheNet prior models, and CellphoneDB input data. Databases
are listed without being loaded, and resources that are not installed or
not configured are reported with their status instead of failing.

## Usage

``` r
ListCCCDB(
  db = c("LIANA", "CellChat", "Nichenetr", "CellphoneDB"),
  species = NULL
)
```

## Arguments

- db:

  CCC methods or resources to inspect. Any subset of \`"LIANA"\`,
  \`"CellChat"\`, \`"Nichenetr"\`, and \`"CellphoneDB"\`; the default
  inspects all supported resources.

- species:

  Optional species filter: any subset of \`"human"\`, \`"mouse"\`, and
  \`"zebrafish"\`.

## Value

A data frame with columns \`db\`, \`resource\`, \`species\`,
\`default\`, \`status\`, and \`description\`. A row with an \`NA\`
resource and a non-\`"available"\` status describes a resource that
could not be inspected.
