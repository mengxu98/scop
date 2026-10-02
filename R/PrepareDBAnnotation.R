#' @title Prepare GO databases
#' @description Prepare Gene Ontology annotations and semantic-similarity data, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character vector selecting `"GO"`, `"GO_BP"`, `"GO_CC"` or `"GO_MF"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`; GO entries may include `semData`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareGO(species = "Homo_sapiens", db = "GO_BP")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareGO <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = c("GO", "GO_BP", "GO_CC", "GO_MF"),
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
    log_message("Unsupported {.arg db} selector for PrepareGO", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      go_categories <- c("GO", "GO_BP", "GO_CC", "GO_MF")
      if (any(db %in% go_categories) &&
        any(!intersect(db, go_categories) %in% names(db_list[[sps]]))
      ) {
        terms <- db[db %in% go_categories]
        bg <- suppressMessages(
          AnnotationDbi::select(
            orgdb,
            keys = AnnotationDbi::keys(orgdb),
            columns = c("GOALL", org_key)
          )
        )
        bg <- unique(bg[
          !is.na(bg[["GOALL"]]),
          c("GOALL", "ONTOLOGYALL", org_key),
          drop = FALSE
        ])
        go_db <- get_namespace_fun("GO.db", "GO.db")
        bg2 <- suppressMessages(
          AnnotationDbi::select(
            go_db,
            keys = AnnotationDbi::keys(go_db),
            columns = c("GOID", "TERM")
          )
        )
        bg <- merge(
          x = bg,
          by.x = "GOALL",
          y = bg2,
          by.y = "GOID",
          all.x = TRUE
        )
        for (subterm in terms) {
          log_message("Preparing database: {.pkg {subterm}}", verbose = verbose)
          if (subterm == "GO") {
            TERM2GENE <- bg[, c("GOALL", org_key)]
            TERM2NAME <- bg[, c("GOALL", "TERM")]
            colnames(TERM2GENE) <- c("Term", default_id_types[[subterm]])
            colnames(TERM2NAME) <- c("Term", "Name")
            TERM2NAME[["ONTOLOGY"]] <- bg[["ONTOLOGYALL"]]
            semData <- NULL
          } else {
            simpleterm <- unlist(strsplit(subterm, split = "_"))[2]
            TERM2GENE <- bg[
              which(bg[["ONTOLOGYALL"]] %in% simpleterm),
              c("GOALL", org_key)
            ]
            TERM2NAME <- bg[
              which(bg[["ONTOLOGYALL"]] %in% simpleterm),
              c("GOALL", "TERM")
            ]
            colnames(TERM2GENE) <- c("Term", default_id_types[[subterm]])
            colnames(TERM2NAME) <- c("Term", "Name")
            TERM2NAME[["ONTOLOGY"]] <- simpleterm
            godata_args <- list(ont = simpleterm)
            if ("annoDb" %in% names(formals(GOSemSim::godata))) {
              godata_args[["annoDb"]] <- orgdb
            } else {
              godata_args[["OrgDb"]] <- orgdb
            }
            semData <- suppressMessages(
              do.call(GOSemSim::godata, godata_args)
            )
          }
          TERM2GENE <- stats::na.omit(unique(TERM2GENE))
          TERM2NAME <- stats::na.omit(unique(TERM2NAME))
          version <- utils::packageVersion(org_sp)
          db_list[[db_species[subterm]]][[subterm]][[
            "TERM2GENE"
          ]] <- TERM2GENE
          db_list[[db_species[subterm]]][[subterm]][[
            "TERM2NAME"
          ]] <- TERM2NAME
          db_list[[db_species[subterm]]][[subterm]][["semData"]] <- semData
          db_list[[db_species[subterm]]][[subterm]][["version"]] <- version
          if (sps == db_species[subterm]) {
            preparedb_cache_annotation(
              db_list[[db_species[subterm]]][[subterm]],
              species = as.character(db_species[subterm]), db = subterm
            )
          }
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare MP databases
#' @description Prepare MP annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"MP"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareMP(species = "Homo_sapiens", db = "MP")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareMP <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "MP",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "MP")) {
    log_message("Unsupported {.arg db} selector for PrepareMP", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (any(db == "MP") && (!"MP" %in% names(db_list[[sps]]))) {
        if (sps != "Mus_musculus") {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the mouse annotation to create the MP database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["MP"] <- "Mus_musculus"
          } else {
            log_message(
              "{.pkg MP} database only support {.val Mus_musculus}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg MP} database", verbose = verbose)
        temp <- tempfile()
        version <- as.character(Sys.Date())
        download(
          quiet = TRUE, url = "https://www.informatics.jax.org/downloads/reports/VOC_MammalianPhenotype.rpt",
          destfile = temp
        )
        mp_name <- utils::read.table(
          temp,
          header = FALSE,
          sep = "\t",
          fill = TRUE,
          quote = ""
        )
        rownames(mp_name) <- mp_name[, 1]
        download(
          quiet = TRUE, url = "https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt",
          destfile = temp
        )
        gene_id <- utils::read.table(
          temp,
          header = FALSE,
          row.names = NULL,
          sep = "\t",
          fill = TRUE,
          quote = ""
        )
        gene_id <- gene_id[, 1:15]
        colnames(gene_id) <- gene_id[1, ]
        gene_id <- gene_id[
          gene_id[, 2] %in% c("Gene", "Pseudogene"), ,
          drop = FALSE
        ]
        rownames(gene_id) <- gene_id[, 1]

        download(
          quiet = TRUE, url = "https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt",
          destfile = temp
        )
        mp_gene <- utils::read.table(
          temp,
          header = FALSE,
          sep = "\t",
          fill = TRUE,
          quote = ""
        )
        mp_gene[["symbol"]] <- gene_id[mp_gene[["V7"]], "3. marker symbol"]
        mp_gene[["MP"]] <- mp_name[mp_gene[, "V5"], 2]
        TERM2GENE <- mp_gene[, c("V5", "symbol")]
        TERM2NAME <- mp_gene[, c("V5", "MP")]

        colnames(TERM2GENE) <- c("Term", default_id_types[["MP"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["MP"]]][["MP"]][["TERM2GENE"]] <- TERM2GENE
        db_list[[db_species["MP"]]][["MP"]][["TERM2NAME"]] <- TERM2NAME
        db_list[[db_species["MP"]]][["MP"]][["version"]] <- version
        if (sps == db_species["MP"]) {
          preparedb_cache_annotation(
            db_list[[db_species["MP"]]][["MP"]],
            species = as.character(db_species["MP"]), db = "MP"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare DO databases
#' @description Prepare DO annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"DO"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareDO(species = "Homo_sapiens", db = "DO")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareDO <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "DO",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "DO")) {
    log_message("Unsupported {.arg db} selector for PrepareDO", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (any(db == "DO") && (!"DO" %in% names(db_list[[sps]]))) {
        log_message("Preparing {.pkg DO} database", verbose = verbose)
        temp <- tempfile(fileext = ".tsv.gz")
        download(
          quiet = TRUE, url = "https://fms.alliancegenome.org/download/DISEASE-ALLIANCE_COMBINED.tsv.gz",
          destfile = temp
        )
        R.utils::gunzip(temp)
        do_all <- utils::read.table(
          gsub(".gz", "", temp),
          header = TRUE,
          sep = "\t",
          fill = TRUE,
          quote = ""
        )
        version <- gsub(
          pattern = ".*Alliance Database Version: ",
          replacement = "",
          x = grep(
            "Alliance Database Version",
            readLines(gsub(".gz", "", temp), warn = FALSE),
            perl = TRUE,
            value = TRUE
          )
        )
        unlink(temp)
        do_sp <- gsub(pattern = "_", replacement = " ", x = sps)
        do_df <- do_all[
          do_all[["DBobjectType"]] == "gene" &
            do_all[["SpeciesName"]] == do_sp, ,
          drop = FALSE
        ]
        if (nrow(do_df) == 0) {
          if (isTRUE(convert_species) && db_species["DO"] != "Homo_sapiens") {
            log_message(
              "Use the human annotation to create the DO database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["DO"] <- "Homo_sapiens"
            do_sp <- gsub(
              pattern = "_",
              replacement = " ",
              x = "Homo_sapiens"
            )
            do_df <- do_all[
              do_all[["DBobjectType"]] == "gene" &
                do_all[["SpeciesName"]] == do_sp, ,
              drop = FALSE
            ]
          } else {
            log_message(
              "Stop the preparation",
              message_type = "error"
            )
          }
        }
        TERM2GENE <- do_df[, c("DOID", "DBObjectSymbol")]
        TERM2NAME <- do_df[, c("DOID", "DOtermName")]
        colnames(TERM2GENE) <- c("Term", default_id_types[["DO"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["DO"]]][["DO"]][["TERM2GENE"]] <- TERM2GENE
        db_list[[db_species["DO"]]][["DO"]][["TERM2NAME"]] <- TERM2NAME
        db_list[[db_species["DO"]]][["DO"]][["version"]] <- version
        if (sps == db_species["DO"]) {
          preparedb_cache_annotation(
            db_list[[db_species["DO"]]][["DO"]],
            species = as.character(db_species["DO"]), db = "DO"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare HPO databases
#' @description Prepare HPO annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"HPO"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareHPO(species = "Homo_sapiens", db = "HPO")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareHPO <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "HPO",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "HPO")) {
    log_message("Unsupported {.arg db} selector for PrepareHPO", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (any(db == "HPO") && (!"HPO" %in% names(db_list[[sps]]))) {
        log_message("Preparing {.pkg HPO} database", verbose = verbose)
        if (!sps %in% c("Homo_sapiens")) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the HPO database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["HPO"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg HPO} database only support {.val Homo_sapiens}. Consider using {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        temp <- tempfile()
        download(
          quiet = TRUE, url = "https://api.github.com/repos/obophenotype/human-phenotype-ontology/releases?per_page=1",
          destfile = temp
        )
        release <- readLines(temp, warn = FALSE)
        release_tag <- regmatches(
          release,
          m = regexpr(
            "(?<=tag_name\\\":\\\")\\S+(?=\\\",\\\"target_commitish)",
            release,
            perl = T
          )
        )
        version <- if (length(release_tag) > 0) {
          release_tag
        } else {
          paste0("Retrieved ", Sys.Date())
        }

        download(
          quiet = TRUE, url = "http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt",
          destfile = temp
        )
        hpo <- utils::read.table(
          temp,
          header = TRUE,
          sep = "\t",
          fill = TRUE,
          quote = ""
        )
        unlink(temp)

        TERM2GENE <- hpo[, c("hpo_id", "gene_symbol")]
        TERM2NAME <- hpo[, c("hpo_id", "hpo_name")]
        colnames(TERM2GENE) <- c("Term", default_id_types[["HPO"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["HPO"]]][["HPO"]][["TERM2GENE"]] <- TERM2GENE
        db_list[[db_species["HPO"]]][["HPO"]][["TERM2NAME"]] <- TERM2NAME
        db_list[[db_species["HPO"]]][["HPO"]][["version"]] <- version
        if (sps == db_species["HPO"]) {
          preparedb_cache_annotation(
            db_list[[db_species["HPO"]]][["HPO"]],
            species = as.character(db_species["HPO"]), db = "HPO"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare PFAM databases
#' @description Prepare PFAM annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"PFAM"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PreparePFAM(species = "Homo_sapiens", db = "PFAM")
#' names(databases[["Homo_sapiens"]])
#' }
PreparePFAM <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "PFAM",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "PFAM")) {
    log_message("Unsupported {.arg db} selector for PreparePFAM", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (any(db == "PFAM") && (!"PFAM" %in% names(db_list[[sps]]))) {
        log_message("Preparing {.pkg PFAM} database", verbose = verbose)
        if (!"PFAM" %in% AnnotationDbi::columns(orgdb)) {
          log_message(
            "{.pkg PFAM} is not in the orgdb: {.val {orgdb}}. Skip this preparation",
            message_type = "warning",
            verbose = verbose
          )
        } else {
          bg <- suppressMessages(
            AnnotationDbi::select(
              orgdb,
              keys = AnnotationDbi::keys(orgdb),
              columns = c("PFAM", org_key)
            )
          )
          bg <- unique(bg[!is.na(bg$PFAM), c("PFAM", org_key), drop = FALSE])
          pfam_de2ac <- get_namespace_fun("PFAM.db", "PFAMDE2AC")
          bg2 <- as.data.frame(
            pfam_de2ac[AnnotationDbi::mappedkeys(pfam_de2ac)]
          )
          rownames(bg2) <- bg2[["ac"]]
          bg[["PFAM_name"]] <- bg2[bg$PFAM, "de"]
          bg[is.na(bg[["PFAM_name"]]), "PFAM_name"] <- bg[
            is.na(bg[["PFAM_name"]]),
            "PFAM"
          ]
          TERM2GENE <- bg[, c("PFAM", org_key)]
          TERM2NAME <- bg[, c("PFAM", "PFAM_name")]
          colnames(TERM2GENE) <- c("Term", default_id_types[["PFAM"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- stats::na.omit(unique(TERM2GENE))
          TERM2NAME <- stats::na.omit(unique(TERM2NAME))
          version <- utils::packageVersion(org_sp)
          db_list[[db_species["PFAM"]]][["PFAM"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["PFAM"]]][["PFAM"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["PFAM"]]][["PFAM"]][["version"]] <- version
          if (sps == db_species["PFAM"]) {
            preparedb_cache_annotation(
              db_list[[db_species["PFAM"]]][["PFAM"]],
              species = as.character(db_species["PFAM"]), db = "PFAM"
            )
          }
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare Chromosome databases
#' @description Prepare Chromosome annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"Chromosome"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareChromosome(species = "Homo_sapiens", db = "Chromosome")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareChromosome <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "Chromosome",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "Chromosome")) {
    log_message("Unsupported {.arg db} selector for PrepareChromosome", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (
        any(db == "Chromosome") && (!"Chromosome" %in% names(db_list[[sps]]))
      ) {
        log_message("Preparing {.pkg Chromosome} database", verbose = verbose)
        orgdbCHR <- get_namespace_fun(org_sp, sub("\\.db$", "CHR", org_sp))
        chr <- as.data.frame(
          orgdbCHR[AnnotationDbi::mappedkeys(orgdbCHR)]
        )
        chr[, 2] <- paste0("chr", chr[, 2])
        TERM2GENE <- chr[, c(2, 1)]
        TERM2NAME <- chr[, c(2, 2)]
        colnames(TERM2GENE) <- c("Term", default_id_types[["Chromosome"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        version <- utils::packageVersion(org_sp)
        db_list[[db_species["Chromosome"]]][["Chromosome"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["Chromosome"]]][["Chromosome"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["Chromosome"]]][["Chromosome"]][[
          "version"
        ]] <- version
        if (sps == db_species["Chromosome"]) {
          preparedb_cache_annotation(
            db_list[[db_species["Chromosome"]]][["Chromosome"]],
            species = as.character(db_species["Chromosome"]), db = "Chromosome"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare GeneType databases
#' @description Prepare GeneType annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"GeneType"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareGeneType(species = "Homo_sapiens", db = "GeneType")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareGeneType <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "GeneType",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "GeneType")) {
    log_message("Unsupported {.arg db} selector for PrepareGeneType", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (any(db == "GeneType") && (!"GeneType" %in% names(db_list[[sps]]))) {
        log_message("Preparing {.pkg GeneType} database", verbose = verbose)
        if (!"GENETYPE" %in% AnnotationDbi::columns(orgdb)) {
          log_message(
            "GENETYPE is not in the orgdb: {.val {org_sp}}. Skip this preparation",
            message_type = "warning",
            verbose = verbose
          )
        } else {
          bg <- suppressMessages(
            AnnotationDbi::select(
              orgdb,
              keys = AnnotationDbi::keys(orgdb),
              columns = c("GENETYPE", org_key)
            )
          )
          TERM2GENE <- bg[, c("GENETYPE", org_key)]
          TERM2NAME <- bg[, c("GENETYPE", "GENETYPE")]
          colnames(TERM2GENE) <- c("Term", default_id_types[["GeneType"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- stats::na.omit(unique(TERM2GENE))
          TERM2NAME <- stats::na.omit(unique(TERM2NAME))
          version <- utils::packageVersion(org_sp)
          db_list[[db_species["GeneType"]]][["GeneType"]][[
            "TERM2GENE"
          ]] <- TERM2GENE
          db_list[[db_species["GeneType"]]][["GeneType"]][[
            "TERM2NAME"
          ]] <- TERM2NAME
          db_list[[db_species["GeneType"]]][["GeneType"]][[
            "version"
          ]] <- version
          if (sps == db_species["GeneType"]) {
            preparedb_cache_annotation(
              db_list[[db_species["GeneType"]]][["GeneType"]],
              species = as.character(db_species["GeneType"]), db = "GeneType"
            )
          }
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare Enzyme databases
#' @description Prepare Enzyme annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"Enzyme"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareEnzyme(species = "Homo_sapiens", db = "Enzyme")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareEnzyme <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "Enzyme",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "Enzyme")) {
    log_message("Unsupported {.arg db} selector for PrepareEnzyme", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (any(db == "Enzyme") && (!"Enzyme" %in% names(db_list[[sps]]))) {
        log_message("Preparing {.pkg Enzyme} database", verbose = verbose)
        if (!"ENZYME" %in% AnnotationDbi::columns(orgdb)) {
          log_message(
            "ENZYME is not in the orgdb: {.val {orgdb}}. Skip this preparation",
            message_type = "warning",
            verbose = verbose
          )
        } else {
          bg <- suppressMessages(
            AnnotationDbi::select(
              orgdb,
              keys = AnnotationDbi::keys(orgdb),
              columns = c("ENZYME", org_key)
            )
          )
          bg1 <- bg2 <- stats::na.omit(bg)
          bg1[, "ENZYME"] <- sapply(
            strsplit(bg1[, "ENZYME"], "\\."),
            function(x) paste0(utils::head(x, 1), collapse = ".")
          )
          bg2[, "ENZYME"] <- sapply(
            strsplit(bg2[, "ENZYME"], "\\."),
            function(x) paste0(utils::head(x, 2), collapse = ".")
          )
          bg <- unique(rbind(bg1, bg2))
          bg[, "ENZYME"] <- gsub(pattern = "\\.-$", "", x = bg[, 2])
          bg[, "ENZYME"] <- paste0("ec:", bg[, "ENZYME"])
          temp <- tempfile()
          download(
            quiet = TRUE, url = "https://ftp.expasy.org/databases/enzyme/enzclass.txt",
            destfile = temp
          )
          enzyme <- utils::read.table(
            temp,
            header = FALSE,
            sep = "\t",
            fill = TRUE,
            quote = ""
          )
          enzyme <- enzyme[
            grep("-.-", enzyme[, 1], fixed = TRUE), ,
            drop = FALSE
          ]
          enzyme <- do.call(rbind, strsplit(enzyme[, 1], split = ". -.-  "))
          enzyme[, 1] <- paste0(
            "ec:",
            gsub(pattern = "( )|(. -)", replacement = "", enzyme[, 1])
          )
          enzyme[, 2] <- gsub(
            pattern = "(^ )|(\\.$)",
            replacement = "",
            enzyme[, 2]
          )
          rownames(enzyme) <- enzyme[, 1]
          for (i in seq_len(nrow(enzyme))) {
            if (grepl(".", enzyme[i, 1], fixed = TRUE)) {
              enzyme[i, 2] <- paste0(
                enzyme[strsplit(enzyme[i, 1], ".", fixed = TRUE)[[1]][1], 2],
                "(",
                enzyme[i, 2],
                ")"
              )
            }
          }
          unlink(temp)
          bg[, "Name"] <- enzyme[bg[, "ENZYME"], 2]
          TERM2GENE <- bg[, c("ENZYME", org_key)]
          TERM2NAME <- bg[, c("ENZYME", "Name")]
          colnames(TERM2GENE) <- c("Term", default_id_types[["Enzyme"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- stats::na.omit(unique(TERM2GENE))
          TERM2NAME <- stats::na.omit(unique(TERM2NAME))
          version <- utils::packageVersion(org_sp)
          db_list[[db_species["Enzyme"]]][["Enzyme"]][[
            "TERM2GENE"
          ]] <- TERM2GENE
          db_list[[db_species["Enzyme"]]][["Enzyme"]][[
            "TERM2NAME"
          ]] <- TERM2NAME
          db_list[[db_species["Enzyme"]]][["Enzyme"]][["version"]] <- version
          if (sps == db_species["Enzyme"]) {
            preparedb_cache_annotation(
              db_list[[db_species["Enzyme"]]][["Enzyme"]],
              species = as.character(db_species["Enzyme"]), db = "Enzyme"
            )
          }
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}
