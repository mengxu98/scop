RunPando <- function(input, config) {
  if (is.null(config$resources$genome)) {
    stop("Pando requires resources$genome matching the chromatin coordinates.", call. = FALSE)
  }
  config$method_args <- utils::modifyList(grn_pando_defaults(), config$method_args)
  if (!identical(config$method_args$interaction_term, ":")) {
    stop("Pando GRN layers require interaction_term=':' to identify TF-region terms.", call. = FALSE)
  }
  if (!identical(as.numeric(config$method_args$extend), 1000000)) {
    stop("Pando 1.1.1 does not forward a custom extend value; supply resources$gene_domains for explicit gene domains.", call. = FALSE)
  }
  native <- grn_pando_fit(input, config)
  grn_pando_group(native, input, list(
    mode = "current_session",
    r_version = as.character(getRversion()),
    r_library = .libPaths()
  ))
}

grn_pando_group <- function(native, input, process = list()) {
  required <- c(
    "coefficients", "goodness_of_fit", "modules", "triplets", "tf_region",
    "region_gene", "gene_sets", "parameters", "package_versions"
  )
  if (!is.list(native) || !all(required %in% names(native))) {
    stop("Pando returned an incomplete native result.", call. = FALSE)
  }
  tables <- list(
    tf_gene = c("modules", "tf", "target", "estimate", "TF", "gene", "aggregated_interaction_coefficient"),
    tf_region = c("tf_region", "tf", "region", "membership", "TF", "region", "membership"),
    region_gene = c("region_gene", "region", "target", "membership", "region", "gene", "membership")
  )
  specs <- lapply(names(tables), function(layer) {
    columns <- tables[[layer]]
    table <- native[[columns[[1L]]]]
    if (!is.data.frame(table) || !all(columns[2:4] %in% names(table))) {
      stop("Pando returned an invalid ", layer, " table.", call. = FALSE)
    }
    membership <- identical(columns[[7L]], "membership")
    list(
      status = "computed", table_path = columns[[1L]],
      source = columns[[2L]], target = columns[[3L]], score = columns[[4L]],
      source_type = columns[[5L]], target_type = columns[[6L]],
      score_type = columns[[7L]], directed = TRUE,
      sign_column = if (membership) NULL else "estimate",
      sign_basis = if (membership) "unknown" else "aggregated_interaction_coefficient",
      rank_by = if (membership) "score" else "absolute_score"
    )
  })
  names(specs) <- names(tables)
  list(
    cells = input$cells,
    native = native,
    layer_specs = specs,
    chains = list(status = "computed", table_path = "triplets"),
    gene_sets = list(status = "computed", table_path = "gene_sets"),
    parameters = native$parameters,
    implementation = "Pando::infer_grn/glm",
    package_versions = native$package_versions,
    execution = utils::modifyList(process, list(status = "completed", artifacts = list()))
  )
}

grn_pando_fit <- function(input, config) {
  args <- config$method_args
  pando <- grn_pando_functions(config$cores)
  genome <- config$resources$genome
  if (is.character(genome) && length(genome) == 1L) {
    if (!requireNamespace(genome, quietly = TRUE)) {
      stop("The selected Pando genome package is not installed: ", genome, call. = FALSE)
    }
    genome <- getExportedValue(genome, genome)
  }
  genome_label <- unique(as.character(GenomeInfoDb::genome(genome)))
  human_hg38 <- any(grepl("hg38|GRCh38", genome_label))
  if (identical(args$peak_to_gene_method, "GREAT") && is.null(config$resources$gene_domains) && !human_hg38) {
    stop("Pando's native GREAT domains use hg38; select peak_to_gene_method='Signac' or supply resources$gene_domains for another genome.", call. = FALSE)
  }
  object <- input$object[, input$cells]
  rna_assay <- config$assay
  peak_assay <- config$atac_assay
  if (!methods::is(object[[peak_assay]], "ChromatinAssay")) {
    stop("Pando requires a ChromatinAssay with gene annotations.", call. = FALSE)
  }
  if (!is.null(config$resources$annotation)) Signac::Annotation(object[[peak_assay]]) <- config$resources$annotation
  if (is.null(Signac::Annotation(object[[peak_assay]]))) {
    stop("Pando requires resources$annotation or a ChromatinAssay gene annotation.", call. = FALSE)
  }
  object <- SeuratObject::SetAssayData(object, assay = rna_assay, layer = "data", new.data = input$expression)
  peak_layer <- config$atac_layer
  peak_expression <- SeuratObject::LayerData(object, assay = peak_assay, layer = peak_layer)
  object <- SeuratObject::SetAssayData(object, assay = peak_assay, layer = "data", new.data = peak_expression)
  SeuratObject::DefaultAssay(object) <- rna_assay
  datasets <- new.env(parent = emptyenv())
  utils::data("motifs", package = "Pando", envir = datasets)
  utils::data("motif2tf", package = "Pando", envir = datasets)
  pfm <- config$resources$motif
  if (is.null(pfm)) pfm <- datasets$motifs
  motif_tfs <- config$resources$motif_tfs
  if (is.null(motif_tfs)) motif_tfs <- datasets$motif2tf
  if (!is.null(input$regulators)) {
    motif_tfs <- motif_tfs[as.character(motif_tfs[[2L]]) %in% input$regulators, , drop = FALSE]
    if (!nrow(motif_tfs)) {
      stop("No Pando motif TFs match the selected regulators.", call. = FALSE)
    }
  }
  regions <- config$resources$regions
  if (is.null(regions)) {
    utils::data("phastConsElements20Mammals.UCSC.hg38", package = "Pando", envir = datasets)
    regions <- datasets$phastConsElements20Mammals.UCSC.hg38
    if (!human_hg38) {
      stop("Pando's bundled conserved regions require hg38; supply resources$regions for another genome.", call. = FALSE)
    }
  }
  template <- get_namespace_fun("Pando", "initiate_grn")(
    object,
    regions = regions,
    peak_assay = peak_assay,
    rna_assay = rna_assay,
    exclude_exons = args$exclude_exons
  )
  template <- get_namespace_fun("Pando", "find_motifs")(
    template,
    pfm = pfm,
    motif_tfs = motif_tfs,
    genome = genome,
    verbose = config$verbose
  )
  motif_ids <- colnames(template@grn@regions@motifs@data)
  tf_map <- get_namespace_fun("Pando", "NetworkTFs")(template)
  if (is.null(dim(tf_map))) {
    tf_names <- intersect(rownames(object[[rna_assay]]), as.character(motif_tfs[[2L]]))
    map_ids <- unique(as.character(motif_tfs[[1L]]))
    tf_map <- Matrix::Matrix(tf_map,
      nrow = length(map_ids), ncol = length(tf_names),
      sparse = TRUE, dimnames = list(map_ids, tf_names)
    )
  }
  if (!identical(rownames(tf_map), motif_ids)) {
    index <- match(motif_ids, rownames(tf_map))
    if (all(is.na(index))) {
      stop("Pando motif matrix and TF map have no matching motif IDs.", call. = FALSE)
    }
    zero <- Matrix::Matrix(0, nrow = 1L, ncol = ncol(tf_map), sparse = TRUE)
    index[is.na(index)] <- nrow(tf_map) + 1L
    tf_map <- rbind(tf_map, zero)[index, , drop = FALSE]
    rownames(tf_map) <- motif_ids
  }
  template@grn@regions@tfs <- tf_map
  genes <- input$targets
  if (is.null(genes)) {
    template@data <- Seurat::FindVariableFeatures(
      template@data,
      assay = rna_assay,
      selection.method = "vst",
      nfeatures = args$nfeatures,
      verbose = FALSE
    )
    genes <- SeuratObject::VariableFeatures(template@data, assay = rna_assay)
  }
  genes <- intersect(as.character(genes), rownames(object[[rna_assay]]))
  if (!length(genes)) {
    stop("Pando has no target genes in the selected RNA assay.", call. = FALSE)
  }
  infer_args <- args[c(
    "peak_to_gene_method", "upstream", "downstream", "extend", "only_tss",
    "tf_cor", "peak_cor", "family", "interaction_term", "adjust_method", "scale"
  )]
  if (!is.null(config$resources$gene_domains)) {
    infer_args$peak_to_gene_domains <- config$resources$gene_domains
  }
  network_name <- "scop_pando"
  fitted <- do.call(pando$infer_grn, c(list(
    object = template,
    genes = genes,
    network_name = network_name,
    method = "glm",
    parallel = config$cores > 1L,
    verbose = config$verbose
  ), infer_args))
  coefficients <- as.data.frame(pando$coef(get_namespace_fun("Pando", "GetNetwork")(fitted)))
  goodness <- as.data.frame(get_namespace_fun("Pando", "gof")(fitted))
  if (!all(c("tf", "region", "target", "estimate", "pval", "padj") %in% names(coefficients)) ||
    !all(c("target", "rsq", "nvariables") %in% names(goodness))) {
    stop("The native Pando coefficient or fit schema changed.", call. = FALSE)
  }
  fitted <- pando$find_modules(
    fitted,
    network = network_name,
    p_thresh = args$p_thresh,
    rsq_thresh = args$rsq_thresh,
    nvar_thresh = args$nvar_thresh,
    min_genes_per_module = args$min_genes_per_module
  )
  native_modules <- get_namespace_fun("Pando", "NetworkModules")(fitted)
  modules <- as.data.frame(native_modules@meta)
  eligible <- as.character(goodness$target[
    is.finite(goodness$rsq) & goodness$rsq > args$rsq_thresh &
      is.finite(goodness$nvariables) & goodness$nvariables > args$nvar_thresh
  ])
  retained <- coefficients[
    coefficients$target %in% eligible & (is.na(coefficients$padj) | coefficients$padj < args$p_thresh), ,
    drop = FALSE
  ]
  if (anyDuplicated(retained[c("tf", "region", "target")])) {
    stop("Pando retained duplicate interaction triplets.", call. = FALSE)
  }
  tf_region <- unique(retained[c("tf", "region")])
  tf_region$membership <- rep.int(1, nrow(tf_region))
  region_gene <- unique(retained[c("region", "target")])
  region_gene$membership <- rep.int(1, nrow(region_gene))
  gene_sets <- c(
    stats::setNames(native_modules@features$genes_pos, paste0(names(native_modules@features$genes_pos), "(+)", recycle0 = TRUE)),
    stats::setNames(native_modules@features$genes_neg, paste0(names(native_modules@features$genes_neg), "(-)", recycle0 = TRUE))
  )
  packages <- c("Pando", "Seurat", "SeuratObject", "Signac")
  list(
    coefficients = coefficients,
    goodness_of_fit = goodness,
    modules = modules,
    triplets = retained,
    tf_region = tf_region,
    region_gene = region_gene,
    gene_sets = gene_sets,
    module_features = native_modules@features,
    cells = input$cells,
    selected_genes = genes,
    parameters = c(args, list(
      rna_assay = rna_assay,
      peak_assay = peak_assay,
      rna_layer = config$layer,
      peak_layer = peak_layer,
      regulators = input$regulators,
      targets = genes,
      method = "glm",
      cores = config$cores,
      seed = config$seed
    )),
    package_versions = vapply(packages, function(package) as.character(getNamespaceVersion(package)), character(1))
  )
}

grn_pando_defaults <- function() {
  list(
    exclude_exons = TRUE, nfeatures = 2000L,
    peak_to_gene_method = "GREAT", upstream = 100000, downstream = 0,
    extend = 1000000, only_tss = FALSE, tf_cor = 0.1, peak_cor = 0,
    family = "gaussian", interaction_term = ":", adjust_method = "fdr",
    scale = FALSE, p_thresh = 0.05, rsq_thresh = 0.1,
    nvar_thresh = 10L, min_genes_per_module = 5L
  )
}

grn_pando_functions <- function(cores) {
  if (!"Pando" %in% loadedNamespaces()) {
    classes <- list(coef = "Network", print = c("Modules", "Network", "Regions"))
    generics <- list2env(list(coef = stats::coef, print = base::print), parent = emptyenv())
    saved <- Map(function(generic, class) {
      stats::setNames(lapply(class, function(name) {
        utils::getS3method(generic, name, optional = TRUE, envir = generics)
      }), class)
    }, names(classes), classes)
    on.exit({
      for (generic in names(saved)) {
        for (class in names(saved[[generic]])) {
          method <- saved[[generic]][[class]]
          if (!is.null(method)) registerS3method(generic, class, method, envir = generics)
        }
      }
    }, add = TRUE)
  }
  package <- if (length(find.package("Pando", quiet = TRUE))) "Pando" else env_r_packages("pando")[[1L]]
  check_r(package, verbose = FALSE)
  namespace <- asNamespace("Pando")
  if (getNamespaceVersion(namespace) != "1.1.1") {
    stop("This Pando adapter requires Pando 1.1.1. See PrepareEnv(modules = 'pando').", call. = FALSE)
  }

  scope <- new.env(parent = namespace)
  upstream <- function(name, args) {
    fun <- get_namespace_fun("Pando", name)
    if (!all(args %in% names(formals(fun)))) {
      stop("Unsupported signature for Pando's ", name, ". See PrepareEnv(modules = 'pando').", call. = FALSE)
    }
    environment(fun) <- scope
    fun
  }
  scope$GRanges <- GenomicRanges::GRanges
  scope$check_if_available <- function(method) invisible(NULL)
  scope$find_peaks_near_genes <- upstream("find_peaks_near_genes", c("peaks", "genes", "method"))
  bind_data <- function(expr) {
    if (!is.call(expr)) return(expr)
    if (identical(expr[[1L]], quote(utils::data))) {
      expr[["package"]] <- "Pando"
      return(expr)
    }
    as.call(lapply(as.list(expr), bind_data))
  }
  body(scope$find_peaks_near_genes) <- bind_data(body(scope$find_peaks_near_genes))
  scope$fit_grn_models <- upstream("fit_grn_models.GRNData", c("object", "genes", "parallel", "method", "..."))
  scope$coef <- upstream("coef.Network", "object")
  scope$find_modules <- upstream("find_modules.Network", c("object", "p_thresh", "rsq_thresh", "nvar_thresh", "min_genes_per_module"))
  if (cores > 1L) {
    scope$map_par <- function(x, fun, parallel, verbose) {
      thisutils::parallelize_fun(x, fun, cores = cores, verbose = verbose)
    }
  }
  list(
    infer_grn = upstream("infer_grn.GRNData", c("object", "genes", "parallel", "method", "...")),
    find_modules = upstream("find_modules.GRNData", c("object", "network", "p_thresh", "rsq_thresh", "nvar_thresh", "min_genes_per_module")),
    coef = scope$coef
  )
}
