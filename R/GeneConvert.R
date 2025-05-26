#' Gene ID conversion function using biomart
#'
#' This function can convert different gene ID types within one species or between two species using the biomart service.
#'
#' @param geneID A vector of the geneID character.
#' @param geneID_from_IDtype Gene ID type of the input \code{geneID}. e.g. "symbol", "ensembl_id", "entrez_id"
#' @param geneID_to_IDtype Gene ID type(s) to convert to. e.g. "symbol", "ensembl_id", "entrez_id"
#' @param species_from Latin names for animals of the input geneID. e.g. "Homo_sapiens","Mus_musculus"
#' @param species_to Latin names for animals of the output geneID. e.g. "Homo_sapiens","Mus_musculus"
#' @param Ensembl_version Ensembl database version. If NULL, use the current release version.
#' @param biomart The name of the BioMart database that you want to connect to.
#' Possible options include "ensembl", "protists_mart", "fungi_mart", and "plants_mart".
#' @param max_tries The maximum number of attempts to connect with the BioMart service.
#' @param mirror Specify an Ensembl mirror to connect to. The valid options here are 'www', 'uswest', 'useast', 'asia'.
#'
#' @return A list with the following elements:
#'   \itemize{
#'     \item \code{geneID_res:} A data.frame contains the all gene IDs mapped in the database with columns: 'from_IDtype','from_geneID','to_IDtype','to_geneID'.
#'     \item \code{geneID_collapse:} The data.frame contains all the successfully converted gene IDs, and the output gene IDs are collapsed into a list. As a result, the 'from_geneID' column (which is set as the row names) of the data.frame is unique.
#'     \item \code{geneID_expand:} The data.frame contains all the successfully converted gene IDs, and the output gene IDs are expanded.
#'     \item \code{Ensembl_version:} Ensembl database version.
#'     \item \code{Datasets:} Datasets available in the selected BioMart database.
#'     \item \code{Attributes:} Attributes available in the selected BioMart database.
#'     \item \code{geneID_unmapped:} A character vector of gene IDs that are unmapped in the database.
#'   }
#'
#' @export
#'
#' @examples
#' res <- GeneConvert(
#'   geneID = c("CDK1", "MKI67", "TOP2A", "AURKA", "CTCF"),
#'   geneID_from_IDtype = "symbol",
#'   geneID_to_IDtype = "entrez_id",
#'   species_from = "Homo_sapiens",
#'   species_to = "Mus_musculus",
#'   Ensembl_version = 103
#' )
#' str(res)
#'
#' # Convert the human genes to mouse homologs,
#' # and replace the raw counts in a Seurat object.
#' data("pancreas_sub")
#' counts <- pancreas_sub@assays$RNA@counts
#' res <- GeneConvert(
#'   geneID = rownames(counts),
#'   geneID_from_IDtype = "symbol",
#'   geneID_to_IDtype = "symbol",
#'   species_from = "Mus_musculus",
#'   species_to = "Homo_sapiens",
#'   Ensembl_version = 103
#' )
#' # Check the number of input and converted gene IDs
#' input_genes <- length(rownames(counts))
#' db_genes <- length(unique(res$geneID_res$from_geneID))
#' converted_genes_input <- length(unique(res$geneID_collapse$from_geneID))
#' converted_genes_output <- length(unique(res$geneID_expand$symbol))
#' cat("Number of input gene IDs:", input_genes, "\n")
#' cat("Number of gene IDs mapped in the database:", db_genes, "\n")
#' cat(
#'   "Number of input gene IDs that were successfully converted:",
#'   converted_genes_input, "\n"
#' )
#' cat("Number of converted gene IDs:", converted_genes_output, "\n")
#'
#' homologs_counts <- stats::aggregate(
#'   x = counts[res$geneID_expand[, "from_geneID"], ],
#'   by = list(res$geneID_expand[, "symbol"]), FUN = sum
#' )
#' rownames(homologs_counts) <- homologs_counts[, 1]
#' homologs_counts <- methods::as(
#'   Matrix::as.matrix(homologs_counts[, -1]),
#'   "dgCMatrix"
#' )
#' homologs_counts
GeneConvert <- function(
    geneID,
    geneID_from_IDtype = "symbol",
    geneID_to_IDtype = "entrez_id",
    species_from = "Homo_sapiens",
    species_to = NULL,
    Ensembl_version = 103,
    biomart = NULL,
    mirror = NULL,
    max_tries = 5) {
  if (requireNamespace("httr", quietly = TRUE)) {
    httr::set_config(httr::config(
      ssl_verifypeer = FALSE,
      ssl_verifyhost = FALSE
    ))
  }

  if (missing(geneID)) {
    stop("'geneID' must be provided.")
  }
  if (is.null(species_to)) {
    species_to <- species_from
  }
  if (is.null(Ensembl_version)) {
    Ensembl_version <- "current_release"
  }
  species_from_split <- unlist(strsplit(species_from, split = "_"))
  species_to_split <- unlist(strsplit(species_to, split = "_"))
  species_from_simp <- paste0(
    tolower(substring(species_from_split[1], 1, 1)),
    species_from_split[2]
  )
  species_to_simp <- paste0(
    tolower(substring(species_to_split[1], 1, 1)),
    species_to_split[2]
  )
  geneID_from_IDtype <- sapply(geneID_from_IDtype, tolower)
  geneID_to_IDtype <- sapply(unique(geneID_to_IDtype), tolower)

  if ("symbol" %in% geneID_from_IDtype) {
    geneID_from_IDtype <- geneID_from_IDtype[geneID_from_IDtype != "symbol"]
    geneID_from_IDtype <- c(
      "ensembl_symbol",
      "entrez_symbol",
      "uniprot_symbol",
      "wiki_symbol",
      geneID_from_IDtype
    )
  }
  from_IDtype <- sapply(
    geneID_from_IDtype,
    function(x) {
      switch(x,
        "ensembl_symbol" = "external_gene_name",
        "ensembl_id" = "ensembl_gene_id",
        "entrez_symbol" = "entrezgene_accession",
        "entrez_id" = "entrezgene_id",
        "uniprot_symbol" = "uniprot_gn_symbol",
        "wiki_symbol" = "wikigene_name",
        x
      )
    }
  )
  names(from_IDtype) <- geneID_from_IDtype

  to_IDtype <- sapply(
    geneID_to_IDtype,
    function(x) {
      switch(x,
        "symbol" = "external_gene_name",
        "ensembl_symbol" = "external_gene_name",
        "entrez_symbol" = "external_gene_name",
        "ensembl_id" = "ensembl_gene_id",
        "entrez_id" = "entrezgene_id",
        x
      )
    }
  )

  if (
    species_from != species_to &&
      all(geneID_to_IDtype %in% c("symbol", "ensembl_id"))
  ) {
    to_IDtype <- sapply(to_IDtype, function(x) {
      switch(x,
        "external_gene_name" = "associated_gene_name",
        "ensembl_gene_id" = "ensembl_gene"
      )
    })
    to_attr <- paste(species_to_simp, to_IDtype, sep = "_homolog_")
    names(to_attr) <- geneID_to_IDtype
  } else {
    to_attr <- to_IDtype
    names(to_attr) <- geneID_to_IDtype
  }

  if (is.null(biomart)) {
    message("Connect to the Ensembl archives...")
    archives <- try_get(
      expr = {
        biomaRt::listEnsemblArchives()
      },
      max_tries = max_tries,
      error_message = "Get errors when connecting with EnsemblArchives..."
    )
    Ensembl_version <- as.character(Ensembl_version)
    if (Ensembl_version == "current_release") {
      url <- archives[which(archives$current_release == "*"), "url"]
      version <- as.character(archives[
        which(archives$current_release == "*"),
        "version"
      ])
      message(
        "Using the ",
        Ensembl_version,
        "(",
        version,
        ")",
        " version of biomart..."
      )
    } else if (Ensembl_version %in% archives$version) {
      url <- archives[which(archives$version == Ensembl_version), "url"]
      version <- as.character(archives[
        which(archives$version == Ensembl_version),
        "version"
      ])
      message("Using the ", version, " version of biomart...")
    } else {
      stop(
        "Ensembl_version is invalid. Must be one of current_release,",
        paste0(archives$version, collapse = ",")
      )
    }

    message("Connecting to the biomart...")
    mart <- try_get(
      expr = {
        if (!is.null(mirror)) {
          biomaRt::useEnsembl(biomart = "ensembl", mirror = mirror)
        } else {
          biomaRt::useMart(biomart = "ensembl", host = url)
        }
      },
      max_tries = max_tries,
      error_message = "Get errors when connecting with ensembl mart..."
    )
    mart_from <- mart_to <- mart
  } else {
    biomart <- match.arg(
      biomart,
      choices = c(
        "ensembl",
        "protists_mart",
        "fungi_mart",
        "plants_mart"
      ),
      several.ok = TRUE
    )
    url <- stats::setNames(
      object = c(
        "https://ensembl.org",
        "https://protists.ensembl.org",
        "https://fungi.ensembl.org",
        "https://plants.ensembl.org"
      ),
      nm = c("ensembl", "protists_mart", "fungi_mart", "plants_mart")
    )
    message("Connecting to the biomart(", biomart, ")...")
    if (length(biomart) == 1) {
      mart <- try_get(
        expr = {
          biomaRt::useMart(biomart = biomart, host = url[biomart])
        },
        max_tries = max_tries,
        error_message = paste0(
          "Get errors when connecting with ensembl mart(",
          biomart,
          ")"
        )
      )
      mart_from <- mart_to <- mart
    } else {
      stop("Supports conversion within one mart only.")
      mart_from <- try_get(
        expr = {
          biomaRt::useMart(
            biomart = biomart[1],
            host = url[biomart[1]]
          )
        },
        max_tries = max_tries,
        error_message = paste0(
          "Get errors when connecting with ensembl mart(",
          biomart[1],
          ")"
        )
      )
      mart_to <- try_get(
        expr = {
          biomaRt::useMart(
            biomart = biomart[2],
            host = url[biomart[2]]
          )
        },
        max_tries = max_tries,
        error_message = paste0(
          "Get errors when connecting with ensembl mart(",
          biomart[2],
          ")"
        )
      )
    }
  }

  message("Searching the dataset ", species_from_simp, " ...")
  Datasets <- try_get(
    expr = {
      biomaRt::listDatasets(mart_from)
    },
    max_tries = max_tries,
    error_message = paste0(
      "Get errors when connecting with ensembl mart(",
      mart_from@biomart,
      ")"
    )
  )
  dataset <- searchDatasets(
    Datasets,
    pattern = species_from_simp
  )[["dataset"]][
    1
  ]
  if (is.null(dataset)) {
    warning(
      paste0(
        "Can not find the dataset for the species: ",
        species_from,
        " (",
        species_from_simp,
        ")"
      ),
      immediate. = TRUE
    )
    return(list(
      geneID_res = NULL,
      geneID_collapse = NULL,
      geneID_expand = NULL,
      Ensembl_version = version,
      Datasets = Datasets,
      Attributes = NULL
    ))
  }

  message("Connecting to the dataset ", dataset, " ...")
  mart1 <- try_get(
    expr = {
      biomaRt::useDataset(
        dataset = dataset,
        mart = mart_from
      )
    },
    max_tries = max_tries,
    error_message = paste0(
      "Get errors when connecting with Dataset(",
      dataset,
      ")"
    )
  )

  if (
    species_from != species_to &&
      any(!geneID_to_IDtype %in% c("symbol", "ensembl_id"))
  ) {
    message("Searching the dataset ", species_to_simp, " ...")
    Datasets2 <- try_get(
      expr = {
        biomaRt::listDatasets(mart_to)
      },
      max_tries = max_tries,
      error_message = paste0(
        "Get errors when connecting with ensembl mart(",
        mart_to@biomart,
        ")"
      )
    )
    dataset2 <- searchDatasets(
      Datasets2,
      pattern = species_to_simp
    )[[
      "dataset"
    ]][1]
    if (is.null(dataset2)) {
      warning(
        paste0(
          "Can not find the dataset for the species: ",
          species_to,
          " (",
          species_to_simp,
          ")"
        ),
        immediate. = TRUE
      )
      return(
        list(
          geneID_res = NULL,
          geneID_collapse = NULL,
          geneID_expand = NULL,
          Ensembl_version = version,
          Datasets = list(Datasets, Datasets2),
          Attributes = NULL
        )
      )
    }

    message("Connecting to the dataset ", dataset2, " ...")
    mart2 <- try_get(
      expr = {
        biomaRt::useDataset(
          dataset = dataset2,
          mart = mart_to
        )
      },
      max_tries = max_tries,
      error_message = paste0(
        "Get errors when connecting with Dataset(",
        dataset2,
        ")"
      )
    )
  }

  Attributes <- biomaRt::listAttributes(mart1)
  from_IDtype <- from_IDtype[from_IDtype %in% Attributes[, "name"]]
  geneID_res_list <- list()
  total <- length(geneID)
  geneID_res <- NULL
  if (any(!to_attr %in% Attributes$name)) {
    to_attr_drop <- to_attr[!to_attr %in% Attributes$name]
    to_attr <- to_attr[to_attr %in% Attributes$name]
    warning(
      paste0(
        "Can not find the attributes for the species ",
        species_from,
        ": ",
        paste(to_attr_drop, collapse = ", ")
      ),
      immediate. = TRUE
    )
    if (length(to_attr) == 0) {
      warning(
        "No attribute found for the species ",
        species_from,
        ". Please check the 'Attributes' in the result.",
        immediate. = TRUE
      )
      return(
        list(
          geneID_res = NULL,
          geneID_collapse = NULL,
          geneID_expand = NULL,
          Ensembl_version = version,
          Datasets = Datasets,
          Attributes = Attributes
        )
      )
    }
  }

  message("Converting the geneIDs...")
  if (species_from != species_to) {
    for (from_attr in from_IDtype) {
      if (length(geneID) > 0) {
        geneID_res1 <- try_get(
          expr = {
            biomaRt::getBM(
              mart = mart1,
              attributes = c(from_attr, "ensembl_gene_id"),
              filters = from_attr,
              values = list(geneID)
            )
          },
          max_tries = max_tries,
          error_message = "Get errors when retrieving information from the BioMart database"
        )

        geneID_res1 <- geneID_res1[,
          c(from_attr, "ensembl_gene_id"),
          drop = FALSE
        ]
        geneID_res1 <- geneID_res1[
          geneID_res1[, from_attr] %in% geneID, ,
          drop = FALSE
        ]
        if (nrow(geneID_res1) == 0) {
          next
        }
        colnames(geneID_res1) <- c("from_geneID", "ensembl_gene_id_tmp")
        from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
        geneID_res1[, "from_IDtype"] <- from_name

        if (all(geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
          geneID_res2 <- try_get(
            expr = {
              biomaRt::getBM(
                mart = mart1,
                attributes = unique(c("ensembl_gene_id", to_attr)),
                filters = "ensembl_gene_id",
                values = list(geneID_res1[, "ensembl_gene_id_tmp"])
              )
            },
            max_tries = max_tries,
            error_message = "Get errors when retrieving information from the BioMart database"
          )
          geneID_res2 <- geneID_res2[, unique(c("ensembl_gene_id", to_attr))]
          geneID_res2 <- geneID_res2[
            geneID_res2[, "ensembl_gene_id"] %in%
              geneID_res1[, "ensembl_gene_id_tmp"], ,
            drop = FALSE
          ]
          if (nrow(geneID_res2) == 0) {
            next
          }
          geneID_res2[, "ensembl_gene_id_tmp"] <- geneID_res2[
            ,
            "ensembl_gene_id"
          ]
          geneID_res2 <- geneID_res2[, c("ensembl_gene_id_tmp", to_attr)]
          geneID_res2 <- unique(
            reshape2::melt(
              geneID_res2,
              id.vars = "ensembl_gene_id_tmp",
              variable.name = "to_IDtype",
              value.name = "to_geneID"
            )
          )
          geneID_res2$to_IDtype <- stats::setNames(
            names(to_attr),
            nm = to_attr
          )[
            geneID_res2$to_IDtype
          ]
          geneID_res_merge <- merge(
            x = geneID_res1,
            y = geneID_res2,
            by = "ensembl_gene_id_tmp"
          )
        } else {
          homolog_ensembl_gene <- paste(
            species_to_simp,
            "ensembl_gene",
            sep = "_homolog_"
          )
          geneID_res2 <- try_get(
            expr = {
              biomaRt::getBM(
                mart = mart1,
                attributes = c("ensembl_gene_id", homolog_ensembl_gene),
                filters = "ensembl_gene_id",
                values = list(geneID_res1[, "ensembl_gene_id_tmp"])
              )
            },
            max_tries = max_tries,
            error_message = "Get errors when retrieving information from the BioMart database"
          )
          geneID_res2 <- geneID_res2[,
            c("ensembl_gene_id", homolog_ensembl_gene),
            drop = FALSE
          ]
          geneID_res2 <- geneID_res2[
            geneID_res2[, "ensembl_gene_id"] %in%
              geneID_res1[, "ensembl_gene_id_tmp"], ,
            drop = FALSE
          ]
          if (nrow(geneID_res2) == 0) {
            next
          }
          colnames(geneID_res2) <- c(
            "ensembl_gene_id_tmp",
            homolog_ensembl_gene
          )

          geneID_res3 <- try_get(
            expr = {
              biomaRt::getBM(
                mart = mart2,
                attributes = unique(c("ensembl_gene_id", to_attr)),
                filters = "ensembl_gene_id",
                values = list(geneID_res2[, homolog_ensembl_gene])
              )
            },
            max_tries = max_tries,
            error_message = "Get errors when retrieving information from the BioMart database"
          )
          geneID_res3 <- geneID_res3[,
            unique(c("ensembl_gene_id", to_attr)),
            drop = FALSE
          ]
          geneID_res3 <- geneID_res3[
            geneID_res3[, "ensembl_gene_id"] %in%
              geneID_res2[, homolog_ensembl_gene], ,
            drop = FALSE
          ]
          if (nrow(geneID_res3) == 0) {
            next
          }
          geneID_res3[, "ensembl_gene_id_tmp2"] <- geneID_res3[
            ,
            "ensembl_gene_id"
          ]
          geneID_res3 <- geneID_res3[, c("ensembl_gene_id_tmp2", to_attr)]
          geneID_res3 <- unique(
            reshape2::melt(
              geneID_res3,
              id.vars = "ensembl_gene_id_tmp2",
              variable.name = "to_IDtype",
              value.name = "to_geneID"
            )
          )
          colnames(geneID_res3)[1] <- homolog_ensembl_gene
          geneID_res3$to_IDtype <- stats::setNames(
            names(to_attr),
            nm = to_attr
          )[
            geneID_res3$to_IDtype
          ]
          geneID_res_merge <- Reduce(
            merge,
            list(geneID_res1, geneID_res2, geneID_res3)
          )
        }
        geneID_res_merge[geneID_res_merge == ""] <- NA
        geneID_res_list[[from_attr]] <- geneID_res_merge[, c(
          "from_IDtype",
          "from_geneID",
          "to_IDtype",
          "to_geneID"
        )]
        ismap <- geneID %in% geneID_res_list[[from_attr]][, "from_geneID"]
        message(paste(sum(ismap), "genes mapped with", from_name))
        geneID <- geneID[!ismap]
      }
    }
  } else {
    for (from_attr in from_IDtype) {
      if (length(geneID) > 0) {
        geneID_res1 <- try_get(
          expr = {
            biomaRt::getBM(
              mart = mart1,
              attributes = unique(
                c("ensembl_gene_id", from_attr, to_attr)
              ),
              filters = from_attr,
              values = list(geneID)
            )
          },
          max_tries = max_tries,
          error_message = "Get errors when retrieving information from the BioMart database"
        )
        geneID_res1 <- geneID_res1[,
          unique(c("ensembl_gene_id", from_attr, to_attr)),
          drop = FALSE
        ]
        geneID_res1 <- geneID_res1[
          geneID_res1[, from_attr] %in% geneID, ,
          drop = FALSE
        ]
        if (nrow(geneID_res1) == 0) {
          next
        }
        geneID_res1[, "ensembl_gene_id_tmp"] <- geneID_res1[, "ensembl_gene_id"]
        geneID_res1_from <- unique(geneID_res1[, c(
          "ensembl_gene_id_tmp",
          from_attr
        )])
        colnames(geneID_res1_from) <- c("ensembl_gene_id_tmp", "from_geneID")
        from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
        geneID_res1_from[, "from_IDtype"] <- from_name

        geneID_res1_to <- unique(geneID_res1[, c(
          "ensembl_gene_id_tmp",
          to_attr
        )])
        geneID_res1_to <- unique(
          reshape2::melt(
            geneID_res1_to,
            id.vars = "ensembl_gene_id_tmp",
            variable.name = "to_IDtype",
            value.name = "to_geneID"
          )
        )
        geneID_res1_to$to_IDtype <- stats::setNames(
          names(to_attr),
          nm = to_attr
        )[
          geneID_res1_to$to_IDtype
        ]

        geneID_res_merge <- merge(
          x = geneID_res1_from,
          y = geneID_res1_to,
          by = "ensembl_gene_id_tmp"
        )

        geneID_res_merge[geneID_res_merge == ""] <- NA
        geneID_res_list[[from_attr]] <- geneID_res_merge[, c(
          "from_IDtype",
          "from_geneID",
          "to_IDtype",
          "to_geneID"
        )]
        ismap <- geneID %in% geneID_res_list[[from_attr]][, "from_geneID"]
        message(paste(sum(ismap), "genes mapped with", from_name))
        geneID <- geneID[!ismap]
      }
    }
  }
  message(
    paste0(
      paste0(rep("=", 30), collapse = ""),
      "\n",
      total - length(geneID),
      " genes mapped\n",
      length(geneID),
      " genes unmapped"
    ),
    "\n",
    paste0(rep("=", 30), collapse = ""),
    "\n"
  )
  geneID_res <- unique(do.call(rbind, geneID_res_list))
  rownames(geneID_res) <- NULL

  if (
    is.null(geneID_res) ||
      nrow(geneID_res) == 0 ||
      all(is.na(geneID_res[["to_geneID"]]))
  ) {
    warning(paste0("None of the gene IDs were converted"), immediate. = TRUE)
    return(list(
      geneID_res = NULL,
      geneID_collapse = NULL,
      geneID_expand = NULL,
      Datasets = Datasets,
      Attributes = Attributes
    ))
  }
  geneID_res_stat <- by(
    geneID_res,
    list(geneID_res[["from_geneID"]], geneID_res[["to_IDtype"]]),
    function(x) {
      data.frame(
        from_IDtype = unique(x[["from_IDtype"]]),
        from_geneID = unique(x[["from_geneID"]]),
        to_IDtype = unique(x[["to_IDtype"]]),
        to_geneID = I(list(unique(x[["to_geneID"]][
          !x[["to_geneID"]] %in% c("", NA)
        ])))
      )
    }
  )
  geneID_collapse <- do.call(rbind, geneID_res_stat)
  geneID_collapse <- unique(as.data.frame(geneID_collapse[, c(
    "from_geneID",
    "to_IDtype",
    "to_geneID"
  )]))
  geneID_collapse <- geneID_collapse[
    sapply(geneID_collapse$to_geneID, length) > 0, ,
    drop = FALSE
  ]
  geneID_collapse <- reshape2::dcast(
    geneID_collapse,
    formula = from_geneID ~ to_IDtype,
    value.var = "to_geneID"
  )
  rownames(geneID_collapse) <- geneID_collapse[, "from_geneID"]
  geneID_expand <- unnest(
    data = geneID_collapse,
    cols = colnames(geneID_collapse)[
      sapply(geneID_collapse, class) %in% c("list", "AsIs")
    ],
    keep_empty = FALSE
  )

  return(
    list(
      geneID_res = geneID_res,
      geneID_collapse = geneID_collapse,
      geneID_expand = geneID_expand,
      Ensembl_version = version,
      Datasets = Datasets,
      Attributes = Attributes,
      geneID_unmapped = geneID
    )
  )
}
