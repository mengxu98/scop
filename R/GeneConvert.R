#' @title Gene ID conversion function using biomart
#'
#' @description
#' This function can convert different gene ID types within one species or between two species using the biomart service.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams PrepareDB
#' @param geneID A vector of the geneID character.
#' @param geneID_from_IDtype Gene ID type of the input `geneID`. e.g. `"symbol"`, `"ensembl_id"`, `"entrez_id"`
#' @param geneID_to_IDtype Gene ID type(s) to convert to. e.g. `"symbol"`, `"ensembl_id"`, `"entrez_id"`.
#' @param species_from Latin names for animals of the input geneID. e.g. `"Homo_sapiens"`, `"Mus_musculus"`.
#' @param species_to Latin names for animals of the output geneID. e.g. `"Homo_sapiens"`, `"Mus_musculus"`.
#' @param biomart The name of the BioMart database that you want to connect to.
#' Possible options include `"ensembl"`, `"protists_mart"`, `"fungi_mart"`, and `"plants_mart"`.
#' @param max_tries The maximum number of attempts to connect with the BioMart service.
#' @param mirror Specify an Ensembl mirror to connect to. The valid options here are `"www"`, `"uswest"`, `"useast"`, `"asia"`.
#'
#' @return A list with the following elements:
#'   \itemize{
#'     \item `geneID_res`: A data.frame contains the all gene IDs mapped in the database with columns: `"from_IDtype"`, `"from_geneID"`, `"to_IDtype"`, `"to_geneID"`.
#'     \item `geneID_collapse`: The data.frame contains all the successfully converted gene IDs, and the output gene IDs are collapsed into a list. As a result, the `"from_geneID"` column (which is set as the row names) of the data.frame is unique.
#'     \item `geneID_expand`: The data.frame contains all the successfully converted gene IDs, and the output gene IDs are expanded.
#'     \item `Ensembl_version`: Ensembl database version.
#'     \item `Datasets`: Datasets available in the selected BioMart database.
#'     \item `Attributes`: Attributes available in the selected BioMart database.
#'     \item `geneID_unmapped`: A character vector of gene IDs that are unmapped in the database.
#'   }
#'
#' @seealso
#' [AnnotateFeatures]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- GeneConvert(
#'   geneID = c("CDK1", "MKI67", "TOP2A", "AURKA", "CTCF"),
#'   species_from = "Homo_sapiens",
#'   species_to = "Mus_musculus"
#' )
#' str(res)
#'
#' # Convert the human genes to mouse homologs,
#' # and replace the raw counts in a Seurat object.
#' data(pancreas_sub)
#' counts <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts"
#' )
#' res <- GeneConvert(
#'   geneID = rownames(counts),
#'   geneID_from_IDtype = "symbol",
#'   geneID_to_IDtype = "symbol",
#'   species_from = "Mus_musculus",
#'   species_to = "Homo_sapiens"
#' )
#'
#' homologs_counts <- stats::aggregate(
#'   x = counts[res$geneID_expand[, "from_geneID"], ],
#'   by = list(res$geneID_expand[, "symbol"]), FUN = sum
#' )
#' rownames(homologs_counts) <- homologs_counts[, 1]
#' homologs_counts <- methods::as(
#'   thisutils::as_matrix(homologs_counts[, -1]),
#'   "dgCMatrix"
#' )
#' homologs_counts[1:5, 1:5]
#' }
GeneConvert <- function(
    geneID,
    geneID_from_IDtype = "symbol",
    geneID_to_IDtype = "entrez_id",
    species_from = "Homo_sapiens",
    species_to = NULL,
    Ensembl_version = NULL,
    biomart = NULL,
    mirror = NULL,
    max_tries = 5,
    verbose = TRUE) {
  if (requireNamespace("httr", quietly = TRUE)) {
    httr::set_config(
      httr::config(
        ssl_verifypeer = FALSE,
        ssl_verifyhost = FALSE
      )
    )
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

  if (species_from != species_to && all(geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
    to_IDtype <- sapply(
      to_IDtype, function(x) {
        switch(x,
          "external_gene_name" = "associated_gene_name",
          "ensembl_gene_id" = "ensembl_gene"
        )
      }
    )
    to_attr <- paste(species_to_simp, to_IDtype, sep = "_homolog_")
    names(to_attr) <- geneID_to_IDtype
  } else {
    to_attr <- to_IDtype
    names(to_attr) <- geneID_to_IDtype
  }

  check_r("biomaRt", verbose = FALSE)
  if (is.null(biomart)) {
    log_message(
      "Connect to the Ensembl archives...",
      verbose = verbose
    )
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
      version <- as.character(
        archives[which(archives$current_release == "*"), "version"]
      )
    } else if (Ensembl_version %in% archives$version) {
      url <- archives[which(archives$version == Ensembl_version), "url"]
      version <- as.character(
        archives[which(archives$version == Ensembl_version), "version"]
      )
    } else {
      log_message(
        "Ensembl_version is invalid. Must be one of {.val {archives$version}}",
        message_type = "error"
      )
    }
    log_message(
      "Using the {.pkg {version}} version of ensembl database...\n",
      "Downloading the ensembl database from {.pkg {url}}...",
      multiline_indent = TRUE,
      verbose = verbose
    )
    mart <- try_get(
      expr = {
        if (!is.null(mirror)) {
          biomaRt::useEnsembl(biomart = "ensembl", mirror = mirror)
        } else {
          biomaRt::useMart(biomart = "ensembl", host = url)
        }
      },
      max_tries = max_tries,
      error_message = "Get errors when connecting with ensembl database..."
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
    log_message(
      "Connecting to the ensembl database...",
      verbose = verbose
    )
    if (length(biomart) == 1) {
      mart <- try_get(
        expr = {
          biomaRt::useMart(biomart = biomart, host = url[biomart])
        },
        max_tries = max_tries,
        error_message = paste0(
          "Get errors when connecting with ensembl database ({.pkg {biomart}})",
        )
      )
      mart_from <- mart_to <- mart
    } else {
      log_message(
        "Supports conversion within one ensembl database only",
        message_type = "error"
      )
      mart_from <- try_get(
        expr = {
          biomaRt::useMart(
            biomart = biomart[1],
            host = url[biomart[1]]
          )
        },
        max_tries = max_tries,
        error_message = paste0(
          "Get errors when connecting with ensembl database ({.pkg {biomart[1]}})",
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
          "Get errors when connecting with ensembl database ({.pkg {biomart[2]}})",
        )
      )
    }
  }

  log_message(
    "Searching the dataset {.pkg {species_from_simp}} ...",
    verbose = verbose
  )
  Datasets <- try_get(
    expr = {
      biomaRt::listDatasets(mart_from)
    },
    max_tries = max_tries,
    error_message = paste0(
      "Get errors when connecting with ensembl database ({.pkg {mart_from@biomart}})",
    )
  )
  dataset <- search_datasets(
    Datasets,
    pattern = species_from_simp,
    verbose = verbose
  )[["dataset"]][1]
  if (is.null(dataset)) {
    log_message(
      "Can not find the dataset for the species: {.pkg {species_from}} ({.pkg {species_from_simp}})",
      message_type = "warning",
      verbose = verbose
    )
    return(
      list(
        geneID_res = NULL,
        geneID_collapse = NULL,
        geneID_expand = NULL,
        Ensembl_version = version,
        Datasets = Datasets,
        Attributes = NULL
      )
    )
  }

  log_message(
    "Connecting to the dataset {.pkg {dataset}} ...",
    verbose = verbose
  )
  mart1 <- try_get(
    expr = {
      biomaRt::useDataset(
        dataset = dataset,
        mart = mart_from
      )
    },
    max_tries = max_tries,
    error_message = paste0(
      "Get errors when connecting with Dataset ({.pkg {dataset}})",
    )
  )

  if (species_from != species_to && any(!geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
    log_message(
      "Searching the dataset {.pkg {species_to_simp}} ...",
      verbose = verbose
    )
    Datasets2 <- try_get(
      expr = {
        biomaRt::listDatasets(mart_to)
      },
      max_tries = max_tries,
      error_message = paste0(
        "Get errors when connecting with ensembl database ({.pkg {mart_to@biomart}})",
      )
    )
    dataset2 <- search_datasets(
      Datasets2,
      pattern = species_to_simp,
      verbose = verbose
    )[["dataset"]][1]
    if (is.null(dataset2)) {
      log_message(
        "Can not find the dataset for the species: {.pkg {species_to}} ({.pkg {species_to_simp}})",
        message_type = "warning",
        verbose = verbose
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

    log_message(
      "Connecting to the dataset {.pkg {dataset2}} ...",
      verbose = verbose
    )
    mart2 <- try_get(
      expr = {
        biomaRt::useDataset(
          dataset = dataset2,
          mart = mart_to
        )
      },
      max_tries = max_tries,
      error_message = paste0(
        "Get errors when connecting with Dataset ({.pkg {dataset2}})",
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
    log_message(
      paste0(
        "Can not find the attributes for the species ",
        species_from,
        ": ",
        paste(to_attr_drop, collapse = ", ")
      ),
      message_type = "warning",
      verbose = verbose
    )
    if (length(to_attr) == 0) {
      log_message(
        "No attribute found for the species {.val {species_from}}. Please check the 'Attributes' in the result",
        message_type = "warning",
        verbose = verbose
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

  log_message("Converting the geneIDs...", verbose = verbose)
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
        log_message(
          "{.val {sum(ismap)}} genes mapped with {.val {from_name}}",
          verbose = verbose
        )
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
        log_message(
          "{.val {sum(ismap)}} genes mapped with {.val {from_name}}",
          verbose = verbose
        )
        geneID <- geneID[!ismap]
      }
    }
  }
  log_message(
    strrep("=", 30), "\n",
    cli::col_green(total - length(geneID), " genes mapped"), "\n",
    cli::col_red(length(geneID), " genes unmapped"), "\n",
    strrep("=", 30),
    verbose = verbose
  )
  geneID_res <- unique(do.call(rbind, geneID_res_list))
  rownames(geneID_res) <- NULL

  if (is.null(geneID_res) || nrow(geneID_res) == 0 || all(is.na(geneID_res[["to_geneID"]]))) {
    log_message(
      "None of the gene IDs were converted",
      message_type = "warning",
      verbose = verbose
    )
    return(
      list(
        geneID_res = NULL,
        geneID_collapse = NULL,
        geneID_expand = NULL,
        Datasets = Datasets,
        Attributes = Attributes
      )
    )
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
  geneID_expand <- unnest_fun(
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

search_datasets <- function(datasets, pattern, verbose = TRUE) {
  col_index <- vapply(
    datasets,
    FUN = function(x) {
      return(
        grepl(
          pattern = pattern,
          x = x,
          ignore.case = TRUE
        )
      )
    },
    FUN.VALUE = logical(length = nrow(datasets))
  )
  row_index <- apply(col_index, 1, any)
  if (any(row_index)) {
    return(datasets[row_index, , drop = FALSE])
  } else {
    log_message(
      "No matching datasets found",
      message_type = "warning",
      verbose = verbose
    )
    return(NULL)
  }
}
