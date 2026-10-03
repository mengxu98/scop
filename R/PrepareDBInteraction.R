#' @title Prepare CORUM databases
#' @description Prepare CORUM annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"CORUM"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareCORUM(species = "Homo_sapiens", db = "CORUM")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareCORUM <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "CORUM",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "CORUM")) {
    log_message("Unsupported {.arg db} selector for PrepareCORUM", message_type = "error")
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
      if (any(db == "CORUM") && (!"CORUM" %in% names(db_list[[sps]]))) {
        if (!sps %in% c("Homo_sapiens")) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the {.pkg CORUM} database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["CORUM"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg CORUM} database only support Homo_sapiens. Consider using convert_species=TRUE",
              message_type = "warning",
              verbose = verbose
            )
            log_message(
              "Stop the preparation",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg CORUM} database", verbose = verbose)
        url <- "https://maayanlab.cloud/static/hdfs/harmonizome/data/corum/gene_set_library_crisp.gmt.gz"
        source_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "CORUM",
          pattern = "^gene_set_library_crisp\\.gmt(\\.gz)?$",
          verbose = verbose
        )
        if (is.null(source_file)) {
          temp <- tempfile(fileext = ".gz")
          download(quiet = TRUE, url = url, destfile = temp)
          R.utils::gunzip(temp)
          source_file <- gsub(".gz", "", temp)
        }
        TERM2GENE <- preparedb_read_gmt_source(source_file)
        version <- "Harmonizome 3.0"
        TERM2NAME <- TERM2GENE[, c(1, 1)]
        colnames(TERM2GENE) <- c("Term", default_id_types[["CORUM"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["CORUM"]]][["CORUM"]][["TERM2GENE"]] <- TERM2GENE
        db_list[[db_species["CORUM"]]][["CORUM"]][["TERM2NAME"]] <- TERM2NAME
        db_list[[db_species["CORUM"]]][["CORUM"]][["version"]] <- version
        if (sps == db_species["CORUM"]) {
          preparedb_cache_annotation(
            db_list[[db_species["CORUM"]]][["CORUM"]],
            species = as.character(db_species["CORUM"]), db = "CORUM"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare TF databases
#' @description Prepare TF annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"TF"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareTF(species = "Homo_sapiens", db = "TF")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareTF <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "TF",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "TF")) {
    log_message("Unsupported {.arg db} selector for PrepareTF", message_type = "error")
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
      if (any(db == "TF") && (!"TF" %in% names(db_list[[sps]]))) {
        log_message("Preparing database: TF", verbose = verbose)

        status <- tryCatch(
          {
            temp <- tempfile()
            url <- paste0(
              "https://raw.githubusercontent.com/mengxu98/datasets/main/AnimalTFDB4/TF_list_final/",
              sps,
              "_TF"
            )
            download(
              quiet = TRUE, url = url,
              destfile = temp
            )
            tf <- utils::read.table(
              temp,
              header = TRUE,
              sep = "\t",
              stringsAsFactors = FALSE,
              fill = TRUE,
              quote = ""
            )
            url <- paste0(
              "https://raw.githubusercontent.com/mengxu98/datasets/main/AnimalTFDB4/Cof_list_final/",
              sps,
              "_Cof"
            )
            download(
              quiet = TRUE, url = url,
              destfile = temp
            )
            tfco <- utils::read.table(
              temp,
              header = TRUE,
              sep = "\t",
              stringsAsFactors = FALSE,
              fill = TRUE,
              quote = ""
            )
            if (!"Symbol" %in% colnames(tf)) {
              if (
                isTRUE(convert_species) && db_species["TF"] != "Homo_sapiens"
              ) {
                log_message(
                  "Use the human annotation to create the TF database for ",
                  sps,
                  message_type = "warning"
                )
                db_species["TF"] <- "Homo_sapiens"
                url <- paste0(
                  "https://raw.githubusercontent.com/mengxu98/datasets/main/AnimalTFDB4/TF_list_final/Homo_sapiens_TF"
                )
                download(quiet = TRUE, url = url, destfile = temp)
                tf <- utils::read.table(
                  temp,
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  fill = TRUE,
                  quote = ""
                )
                url <- paste0(
                  "https://raw.githubusercontent.com/mengxu98/datasets/main/AnimalTFDB4/Cof_list_final/Homo_sapiens_Cof"
                )
                download(quiet = TRUE, url = url, destfile = temp)
                tfco <- utils::read.table(
                  temp,
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  fill = TRUE,
                  quote = ""
                )
              } else {
                log_message(
                  "Stop the preparation.",
                  message_type = "error"
                )
              }
            }
            unlink(temp)
            version <- "AnimalTFDB4"
          },
          error = identity
        )

        if (inherits(status, "error")) {
          temp <- tempfile()
          url <- paste0(
            "https://raw.githubusercontent.com/GuoBioinfoLab/AnimalTFDB3/master/AnimalTFDB3/static/AnimalTFDB3/download/",
            sps,
            "_TF"
          )
          download(quiet = TRUE, url = url, destfile = temp)
          tf <- utils::read.table(
            temp,
            header = TRUE,
            sep = "\t",
            stringsAsFactors = FALSE,
            fill = TRUE,
            quote = ""
          )
          url <- paste0(
            "https://raw.githubusercontent.com/GuoBioinfoLab/AnimalTFDB3/master/AnimalTFDB3/static/AnimalTFDB3/download/",
            sps,
            "_TF_cofactors"
          )
          download(quiet = TRUE, url = url, destfile = temp)
          tfco <- utils::read.table(
            temp,
            header = TRUE,
            sep = "\t",
            stringsAsFactors = FALSE,
            fill = TRUE,
            quote = ""
          )
          if (!"Symbol" %in% colnames(tf)) {
            if (isTRUE(convert_species) && db_species["TF"] != "Homo_sapiens") {
              log_message(
                "Use the human annotation to create the TF database for {.val {sps}}",
                message_type = "warning"
              )
              db_species["TF"] <- "Homo_sapiens"
              url <- c(
                "https://raw.githubusercontent.com/GuoBioinfoLab/AnimalTFDB3/master/AnimalTFDB3/static/AnimalTFDB3/download/Homo_sapiens_TF"
              )
              download(quiet = TRUE, url = url, destfile = temp)
              tf <- utils::read.table(
                temp,
                header = TRUE,
                sep = "\t",
                stringsAsFactors = FALSE,
                fill = TRUE,
                quote = ""
              )
              url <- paste0(
                "https://raw.githubusercontent.com/GuoBioinfoLab/AnimalTFDB3/master/AnimalTFDB3/static/AnimalTFDB3/download/Homo_sapiens_TF_cofactors"
              )
              download(
                quiet = TRUE, url = url,
                destfile = temp
              )
              tfco <- utils::read.table(
                temp,
                header = TRUE,
                sep = "\t",
                stringsAsFactors = FALSE,
                fill = TRUE,
                quote = ""
              )
            } else {
              log_message(
                "Stop the preparation",
                message_type = "error"
              )
            }
          }
          unlink(temp)
          version <- "AnimalTFDB3"
        }

        TERM2GENE <- rbind(
          data.frame("Term" = "TF", "symbol" = tf[["Symbol"]]),
          data.frame("Term" = "TF cofactor", "symbol" = tfco[["Symbol"]])
        )
        TERM2NAME <- data.frame(
          "Term" = c("TF", "TF cofactor"),
          "Name" = c("TF", "TF cofactor")
        )
        colnames(TERM2GENE) <- c("Term", default_id_types[["TF"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["TF"]]][["TF"]][["TERM2GENE"]] <- TERM2GENE
        db_list[[db_species["TF"]]][["TF"]][["TERM2NAME"]] <- TERM2NAME
        db_list[[db_species["TF"]]][["TF"]][["version"]] <- version
        if (sps == db_species["TF"]) {
          preparedb_cache_annotation(
            db_list[[db_species["TF"]]][["TF"]],
            species = as.character(db_species["TF"]), db = "TF"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare CSPA databases
#' @description Prepare CSPA annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"CSPA"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareCSPA(species = "Homo_sapiens", db = "CSPA")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareCSPA <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "CSPA",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "CSPA")) {
    log_message("Unsupported {.arg db} selector for PrepareCSPA", message_type = "error")
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
      if (any(db == "CSPA") && (!"CSPA" %in% names(db_list[[sps]]))) {
        if (!sps %in% c("Homo_sapiens", "Mus_musculus")) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the CSPA database for {.val {sps}}",
              message_type = "warning"
            )
            db_species["CSPA"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg CSPA} database only support {.val {c('Homo_sapiens', 'Mus_musculus')}}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        check_r("openxlsx", verbose = FALSE)
        log_message("Preparing database: CSPA", verbose = verbose)
        url <- "https://raw.githubusercontent.com/mengxu98/datasets/main/CSPA/S1_File.xlsx"
        source_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "CSPA",
          pattern = "^S1_File\\.xlsx$",
          verbose = verbose
        )
        source_is_temp <- FALSE
        if (is.null(source_file)) {
          source_file <- tempfile(fileext = ".xlsx")
          source_is_temp <- TRUE
          download(
            quiet = TRUE, url = url,
            destfile = source_file,
            mode = "wb"
          )
        }
        surfacepro <- get_namespace_fun(
          "openxlsx", "read.xlsx"
        )(source_file, sheet = 1)
        if (isTRUE(source_is_temp)) {
          unlink(source_file)
        }
        surfacepro <- surfacepro[
          surfacepro[["organism"]] ==
            switch(db_species["CSPA"],
              "Homo_sapiens" = "Human",
              "Mus_musculus" = "Mouse"
            ), ,
          drop = FALSE
        ]
        TERM2GENE <- data.frame(
          "Term" = "SurfaceProtein",
          "symbol" = surfacepro[["ENTREZ.gene.symbol"]]
        )
        TERM2NAME <- data.frame(
          "Term" = "SurfaceProtein",
          "Name" = "SurfaceProtein"
        )
        colnames(TERM2GENE) <- c("Term", default_id_types[["CSPA"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        version <- "CSPA"
        db_list[[db_species["CSPA"]]][["CSPA"]][["TERM2GENE"]] <- TERM2GENE
        db_list[[db_species["CSPA"]]][["CSPA"]][["TERM2NAME"]] <- TERM2NAME
        db_list[[db_species["CSPA"]]][["CSPA"]][["version"]] <- version
        if (sps == db_species["CSPA"]) {
          preparedb_cache_annotation(
            db_list[[db_species["CSPA"]]][["CSPA"]],
            species = as.character(db_species["CSPA"]), db = "CSPA"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare Surfaceome databases
#' @description Prepare Surfaceome annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"Surfaceome"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareSurfaceome(species = "Homo_sapiens", db = "Surfaceome")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareSurfaceome <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "Surfaceome",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "Surfaceome")) {
    log_message("Unsupported {.arg db} selector for PrepareSurfaceome", message_type = "error")
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
        any(db == "Surfaceome") && (!"Surfaceome" %in% names(db_list[[sps]]))
      ) {
        if (!sps %in% c("Homo_sapiens")) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the Surfaceome database for {.val {sps}}",
              message_type = "warning"
            )
            db_species["Surfaceome"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg Surfaceome} database only support {.val Homo_sapiens}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        check_r("openxlsx", verbose = FALSE)
        log_message("Preparing database: Surfaceome", verbose = verbose)
        url <- "http://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx"
        source_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "Surfaceome",
          pattern = "^table_S3_surfaceome\\.xlsx$",
          verbose = verbose
        )
        source_is_temp <- FALSE
        if (is.null(source_file)) {
          source_file <- tempfile(fileext = ".xlsx")
          source_is_temp <- TRUE
          download(
            quiet = TRUE, url = url,
            destfile = source_file,
            mode = ifelse(.Platform$OS.type == "windows", "wb", "w")
          )
        }
        surfaceome <- get_namespace_fun("openxlsx", "read.xlsx")(
          source_file,
          sheet = 2,
          colNames = TRUE,
          startRow = 2
        )
        if (isTRUE(source_is_temp)) {
          unlink(source_file)
        }
        TERM2GENE <- data.frame(
          "Term" = "SurfaceProtein",
          "symbol" = surfaceome[["UniProt.gene"]]
        )
        TERM2NAME <- data.frame(
          "Term" = "SurfaceProtein",
          "Name" = "SurfaceProtein"
        )
        colnames(TERM2GENE) <- c("Term", default_id_types[["Surfaceome"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        version <- "Surfaceome"
        db_list[[db_species["Surfaceome"]]][["Surfaceome"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["Surfaceome"]]][["Surfaceome"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["Surfaceome"]]][["Surfaceome"]][[
          "version"
        ]] <- version
        if (sps == db_species["Surfaceome"]) {
          preparedb_cache_annotation(
            db_list[[db_species["Surfaceome"]]][["Surfaceome"]],
            species = as.character(db_species["Surfaceome"]), db = "Surfaceome"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare SPRomeDB databases
#' @description Prepare SPRomeDB annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"SPRomeDB"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareSPRomeDB(species = "Homo_sapiens", db = "SPRomeDB")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareSPRomeDB <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "SPRomeDB",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "SPRomeDB")) {
    log_message("Unsupported {.arg db} selector for PrepareSPRomeDB", message_type = "error")
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
      if (any(db == "SPRomeDB") && (!"SPRomeDB" %in% names(db_list[[sps]]))) {
        if (!sps %in% c("Homo_sapiens")) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the SPRomeDB database for {.val {sps}}",
              message_type = "warning"
            )
            db_species["SPRomeDB"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg SPRomeDB} database only support {.val Homo_sapiens}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg SPRomeDB} database", verbose = verbose)
        url <- "http://119.3.41.228/SPRomeDB/files/download/secreted_proteins_SPRomeDB.csv"
        source_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "SPRomeDB",
          pattern = "^secreted_proteins_SPRomeDB\\.csv$",
          verbose = verbose
        )
        source_is_temp <- FALSE
        if (is.null(source_file)) {
          source_file <- tempfile()
          source_is_temp <- TRUE
          download(quiet = TRUE, url = url, destfile = source_file)
        }
        spromedb <- utils::read.csv(source_file, header = TRUE)
        if (isTRUE(source_is_temp)) {
          unlink(source_file)
        }
        TERM2GENE <- data.frame(
          "Term" = "SecretoryProtein",
          "entrez_id" = unlist(strsplit(spromedb$Gene_ID, ";"))
        )
        TERM2NAME <- data.frame(
          "Term" = "SecretoryProtein",
          "Name" = "SecretoryProtein"
        )
        colnames(TERM2GENE) <- c("Term", default_id_types[["SPRomeDB"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        version <- "SPRomeDB"
        db_list[[db_species["SPRomeDB"]]][["SPRomeDB"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["SPRomeDB"]]][["SPRomeDB"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["SPRomeDB"]]][["SPRomeDB"]][[
          "version"
        ]] <- version
        if (sps == db_species["SPRomeDB"]) {
          preparedb_cache_annotation(
            db_list[[db_species["SPRomeDB"]]][["SPRomeDB"]],
            species = as.character(db_species["SPRomeDB"]), db = "SPRomeDB"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare VerSeDa databases
#' @description Prepare VerSeDa annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"VerSeDa"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareVerSeDa(species = "Homo_sapiens", db = "VerSeDa")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareVerSeDa <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "VerSeDa",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "VerSeDa")) {
    log_message("Unsupported {.arg db} selector for PrepareVerSeDa", message_type = "error")
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
      if (any(db == "VerSeDa") && (!"VerSeDa" %in% names(db_list[[sps]]))) {
        temp <- tempfile()
        download(
          quiet = TRUE, url = "http://genomics.cicbiogune.es/VerSeDa/downloads.php",
          destfile = temp
        )
        verseda_sps <- readLines(temp)
        verseda_sps <- regmatches(
          verseda_sps,
          m = regexpr(
            "(?<=Downloads/)\\S+(?=\\.zip)",
            verseda_sps,
            perl = TRUE
          )
        )
        verseda_sps <- setdiff(
          verseda_sps,
          c("NonRefined", "NonRefined_Curated", "Refined", "Refined_Curated")
        )
        if (!tolower(sps) %in% verseda_sps) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the VerSeDa database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["VerSeDa"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg VerSeDa} database only support {.val {verseda_sps}}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg VerSeDa} database", verbose = verbose)
        temp <- tempfile(fileext = ".zip")
        url <- paste0(
          "http://genomics.cicbiogune.es/VerSeDa/Downloads/",
          tolower(db_species["VerSeDa"]),
          ".zip"
        )
        download(quiet = TRUE, url = url, destfile = temp)
        con <- unz(
          temp,
          paste0(
            tolower(db_species["VerSeDa"]),
            "/",
            tolower(db_species["VerSeDa"]),
            "_Refined.sequences"
          )
        )
        verseda <- readLines(con)
        close(con)
        unlink(temp)
        verseda <- verseda[grep("^>", verseda)]
        verseda <- gsub("^>|\\.\\d+", "", verseda)
        verseda_id <- GeneConvert(
          geneID = verseda,
          geneID_from_IDtype = c(
            "ensembl_peptide_id",
            "refseq_peptide",
            "refseq_peptide_predicted",
            "uniprot_isoform",
            "uniprotswissprot",
            "uniprotsptrembl"
          ),
          geneID_to_IDtype = "symbol",
          species_from = db_species["VerSeDa"]
        )
        TERM2GENE <- data.frame(
          "Term" = "SecretoryProtein",
          "symbol" = unique(verseda_id$geneID_expand$symbol)
        )
        TERM2NAME <- data.frame(
          "Term" = "SecretoryProtein",
          "Name" = "SecretoryProtein"
        )
        colnames(TERM2GENE) <- c("Term", default_id_types[["VerSeDa"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        version <- "VerSeDa"
        db_list[[db_species["VerSeDa"]]][["VerSeDa"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["VerSeDa"]]][["VerSeDa"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["VerSeDa"]]][["VerSeDa"]][["version"]] <- version
        if (sps == db_species["VerSeDa"]) {
          preparedb_cache_annotation(
            db_list[[db_species["VerSeDa"]]][["VerSeDa"]],
            species = as.character(db_species["VerSeDa"]), db = "VerSeDa"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare TFLink databases
#' @description Prepare TFLink annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"TFLink"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareTFLink(species = "Homo_sapiens", db = "TFLink")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareTFLink <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "TFLink",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "TFLink")) {
    log_message("Unsupported {.arg db} selector for PrepareTFLink", message_type = "error")
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
      if (any(db == "TFLink") && (!"TFLink" %in% names(db_list[[sps]]))) {
        tflink_sp <- c(
          "Homo_sapiens",
          "Mus_musculus",
          "Rattus_norvegicus",
          "Danio_rerio",
          "Drosophila_melanogaster",
          "Caenorhabditis_elegans",
          "Saccharomyces_cerevisiae"
        )
        if (!sps %in% tflink_sp) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the TFLink database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["TFLink"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg TFLink} database only support {.val {tflink_sp}}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg TFLink} database", verbose = verbose)
        url <- paste0(
          "https://cdn.netbiol.org/tflink/download_files/TFLink_",
          db_species["TFLink"],
          "_interactions_All_GMT_proteinName_v1.0.gmt"
        )
        source_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "TFLink",
          pattern = paste0(
            "^TFLink_",
            db_species[["TFLink"]],
            "_interactions_All_GMT_proteinName_v1\\.0\\.gmt$"
          ),
          verbose = verbose
        )
        source_is_temp <- FALSE
        if (is.null(source_file)) {
          source_file <- tempfile()
          source_is_temp <- TRUE
          download(quiet = TRUE, url = url, destfile = source_file)
        }
        TERM2GENE <- clusterProfiler::read.gmt(source_file)
        if (isTRUE(source_is_temp)) {
          unlink(source_file)
        }
        version <- "v1.0"
        TERM2NAME <- TERM2GENE[, c(1, 1)]
        colnames(TERM2GENE) <- c("Term", default_id_types[["TFLink"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["TFLink"]]][["TFLink"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["TFLink"]]][["TFLink"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["TFLink"]]][["TFLink"]][["version"]] <- version
        if (sps == db_species["TFLink"]) {
          preparedb_cache_annotation(
            db_list[[db_species["TFLink"]]][["TFLink"]],
            species = as.character(db_species["TFLink"]), db = "TFLink"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare hTFtarget databases
#' @description Prepare hTFtarget annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"hTFtarget"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareHTFtarget(species = "Homo_sapiens", db = "hTFtarget")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareHTFtarget <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "hTFtarget",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "hTFtarget")) {
    log_message("Unsupported {.arg db} selector for PrepareHTFtarget", message_type = "error")
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
        any(db == "hTFtarget") && (!"hTFtarget" %in% names(db_list[[sps]]))
      ) {
        if (!sps %in% "Homo_sapiens") {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the hTFtarget database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["hTFtarget"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg hTFtarget} database only support {.val Homo_sapiens}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg hTFtarget} database", verbose = verbose)
        url <- "https://guolab.wchscu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"
        source_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "hTFtarget",
          pattern = "^tf-target-infomation\\.txt$",
          verbose = verbose
        )
        source_is_temp <- FALSE
        if (is.null(source_file)) {
          source_file <- tempfile()
          source_is_temp <- TRUE
          download(quiet = TRUE, url = url, destfile = source_file)
        }
        TERM2GENE <- utils::read.table(source_file, header = TRUE, fill = T, sep = "\t")
        if (isTRUE(source_is_temp)) {
          unlink(source_file)
        }
        version <- "v1.0"
        TERM2NAME <- TERM2GENE[, c(1, 1)]
        colnames(TERM2GENE) <- c("Term", default_id_types[["hTFtarget"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["hTFtarget"]]][["hTFtarget"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["hTFtarget"]]][["hTFtarget"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["hTFtarget"]]][["hTFtarget"]][[
          "version"
        ]] <- version
        if (sps == db_species["hTFtarget"]) {
          preparedb_cache_annotation(
            db_list[[db_species["hTFtarget"]]][["hTFtarget"]],
            species = as.character(db_species["hTFtarget"]), db = "hTFtarget"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare TRRUST databases
#' @description Prepare TRRUST annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"TRRUST"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareTRRUST(species = "Homo_sapiens", db = "TRRUST")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareTRRUST <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "TRRUST",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "TRRUST")) {
    log_message("Unsupported {.arg db} selector for PrepareTRRUST", message_type = "error")
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
      if (any(db == "TRRUST") && (!"TRRUST" %in% names(db_list[[sps]]))) {
        if (!sps %in% c("Homo_sapiens", "Mus_musculus")) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the TRRUST database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["TRRUST"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg TRRUST} database only support {.val {c('Homo_sapiens', 'Mus_musculus')}}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg TRRUST} database", verbose = verbose)
        url <- switch(db_species["TRRUST"],
          "Homo_sapiens" = "https://raw.githubusercontent.com/bioinfonerd/Transcription-Factor-Databases/master/Ttrust_v2/trrust_rawdata.human.tsv",
          "Mus_musculus" = "https://raw.githubusercontent.com/bioinfonerd/Transcription-Factor-Databases/master/Ttrust_v2/trrust_rawdata.mouse.tsv.gz"
        )
        trrust_file <- switch(db_species["TRRUST"],
          "Homo_sapiens" = "trrust_rawdata.human.tsv",
          "Mus_musculus" = "trrust_rawdata.mouse.tsv.gz"
        )
        source_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "TRRUST",
          pattern = switch(db_species["TRRUST"],
            "Homo_sapiens" = "^trrust_rawdata\\.human\\.tsv$",
            "Mus_musculus" = "^trrust_rawdata\\.mouse\\.tsv\\.gz$"
          ),
          verbose = verbose
        )
        source_is_temp <- FALSE
        if (is.null(source_file)) {
          source_file <- tempfile(fileext = ifelse(endsWith(url, "gz"), ".gz", ""))
          source_is_temp <- TRUE
          download(quiet = TRUE, url = url, destfile = source_file)
        }
        if (grepl("\\.gz$", source_file, ignore.case = TRUE)) {
          temp <- tempfile(fileext = ".gz")
          file.copy(source_file, temp, overwrite = TRUE)
          R.utils::gunzip(temp)
          TERM2GENE <- utils::read.table(
            gsub(".gz$", "", temp),
            header = FALSE,
            fill = T,
            sep = "\t"
          )[, 1:2]
          unlink(gsub(".gz$", "", temp))
        } else {
          TERM2GENE <- utils::read.table(
            source_file,
            header = FALSE,
            fill = T,
            sep = "\t"
          )[, 1:2]
        }
        if (isTRUE(source_is_temp)) {
          unlink(source_file)
        }
        version <- "v2.0"
        TERM2NAME <- TERM2GENE[, c(1, 1)]
        colnames(TERM2GENE) <- c("Term", default_id_types[["TRRUST"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["TRRUST"]]][["TRRUST"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["TRRUST"]]][["TRRUST"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["TRRUST"]]][["TRRUST"]][["version"]] <- version
        if (sps == db_species["TRRUST"]) {
          preparedb_cache_annotation(
            db_list[[db_species["TRRUST"]]][["TRRUST"]],
            species = as.character(db_species["TRRUST"]), db = "TRRUST"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare JASPAR databases
#' @description Prepare JASPAR annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"JASPAR"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareJASPAR(species = "Homo_sapiens", db = "JASPAR")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareJASPAR <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "JASPAR",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "JASPAR")) {
    log_message("Unsupported {.arg db} selector for PrepareJASPAR", message_type = "error")
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
      if (any(db == "JASPAR") && (!"JASPAR" %in% names(db_list[[sps]]))) {
        if (!sps %in% c("Homo_sapiens")) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the JASPAR database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["JASPAR"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg JASPAR} database only support {.val Homo_sapiens}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg JASPAR} database", verbose = verbose)
        url <- "https://maayanlab.cloud/static/hdfs/harmonizome/data/jasparpwm/gene_set_library_crisp.gmt.gz"
        source_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "JASPAR",
          pattern = "^gene_set_library_crisp\\.gmt(\\.gz)?$",
          verbose = verbose
        )
        if (is.null(source_file)) {
          temp <- tempfile(fileext = ".gz")
          download(quiet = TRUE, url = url, destfile = temp)
          R.utils::gunzip(temp)
          source_file <- gsub(".gz", "", temp)
        }
        TERM2GENE <- preparedb_read_gmt_source(source_file)
        version <- "Harmonizome 3.0"
        TERM2NAME <- TERM2GENE[, c(1, 1)]
        colnames(TERM2GENE) <- c("Term", default_id_types[["JASPAR"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["JASPAR"]]][["JASPAR"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["JASPAR"]]][["JASPAR"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["JASPAR"]]][["JASPAR"]][["version"]] <- version
        if (sps == db_species["JASPAR"]) {
          preparedb_cache_annotation(
            db_list[[db_species["JASPAR"]]][["JASPAR"]],
            species = as.character(db_species["JASPAR"]), db = "JASPAR"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare ENCODE databases
#' @description Prepare ENCODE annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"ENCODE"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareENCODE(species = "Homo_sapiens", db = "ENCODE")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareENCODE <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "ENCODE",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "ENCODE")) {
    log_message("Unsupported {.arg db} selector for PrepareENCODE", message_type = "error")
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
      if (any(db == "ENCODE") && (!"ENCODE" %in% names(db_list[[sps]]))) {
        if (!sps %in% c("Homo_sapiens")) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the ENCODE database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["ENCODE"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg ENCODE} database only support {.val Homo_sapiens}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg ENCODE} database", verbose = verbose)
        url <- "https://maayanlab.cloud/static/hdfs/harmonizome/data/encodetfppi/gene_set_library_crisp.gmt.gz"
        source_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "ENCODE",
          pattern = "^gene_set_library_crisp\\.gmt(\\.gz)?$",
          verbose = verbose
        )
        if (is.null(source_file)) {
          temp <- tempfile(fileext = ".gz")
          download(quiet = TRUE, url = url, destfile = temp)
          R.utils::gunzip(temp)
          source_file <- gsub(".gz", "", temp)
        }
        TERM2GENE <- preparedb_read_gmt_source(source_file)
        version <- "Harmonizome 3.0"
        TERM2NAME <- TERM2GENE[, c(1, 1)]
        colnames(TERM2GENE) <- c("Term", default_id_types[["ENCODE"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["ENCODE"]]][["ENCODE"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["ENCODE"]]][["ENCODE"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["ENCODE"]]][["ENCODE"]][["version"]] <- version
        if (sps == db_species["ENCODE"]) {
          preparedb_cache_annotation(
            db_list[[db_species["ENCODE"]]][["ENCODE"]],
            species = as.character(db_species["ENCODE"]), db = "ENCODE"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare CellTalk databases
#' @description Prepare CellTalk annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"CellTalk"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareCellTalk(species = "Homo_sapiens", db = "CellTalk")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareCellTalk <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "CellTalk",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "CellTalk")) {
    log_message("Unsupported {.arg db} selector for PrepareCellTalk", message_type = "error")
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
      ccc_db_use <- intersect(db, c("CellTalk", "CellChat"))
      if (length(ccc_db_use) > 0L &&
        any(!ccc_db_use %in% names(db_list[[sps]]))) {
        ccc_prepared <- PrepareCCCDB(
          species = sps,
          db = ccc_db_use,
          convert_species = convert_species,
          data_dir = data_dir,
          db_version = db_version,
          db_update = db_update,
          verbose = verbose
        )
        for (ccc_db in ccc_db_use) {
          if (!ccc_db %in% names(ccc_prepared[[sps]])) next
          db_list[[sps]][[ccc_db]] <- ccc_prepared[[sps]][[ccc_db]]
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare CellChat databases
#' @description Prepare CellChat annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"CellChat"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareCellChat(species = "Homo_sapiens", db = "CellChat")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareCellChat <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "CellChat",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "CellChat")) {
    log_message("Unsupported {.arg db} selector for PrepareCellChat", message_type = "error")
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
      ccc_db_use <- intersect(db, c("CellTalk", "CellChat"))
      if (length(ccc_db_use) > 0L &&
        any(!ccc_db_use %in% names(db_list[[sps]]))) {
        ccc_prepared <- PrepareCCCDB(
          species = sps,
          db = ccc_db_use,
          convert_species = convert_species,
          data_dir = data_dir,
          db_version = db_version,
          db_update = db_update,
          verbose = verbose
        )
        for (ccc_db in ccc_db_use) {
          if (!ccc_db %in% names(ccc_prepared[[sps]])) next
          db_list[[sps]][[ccc_db]] <- ccc_prepared[[sps]][[ccc_db]]
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}
