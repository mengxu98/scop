#' @title Run scOMM label prediction
#'
#' @description
#' Run `scOMM` on shared features between a reference object and a query object,
#' write predicted labels and class scores into query metadata, and optionally
#' evaluate predictions against a truth label.
#'
#' @md
#' @inheritParams standard_scop
#' @param reference Reference `Seurat` object used for supervision.
#' @param reference_assay Assay used in the reference object.
#' @param query_assay Assay used in the query object.
#' @param reference_label Metadata column in the reference used as supervision labels.
#' @param features Shared features passed to `scOMM`. If `NULL`, reference variable
#' features are used.
#' @param prediction_prefix Prefix added to prediction metadata columns. The
#' default creates `scomm_prediction`, `scomm_score.<class>`, and
#' `scomm_score.max`.
#' @param evaluate Whether to compute prediction metrics against a truth label.
#' @param truth_col Metadata column in `srt` used as the truth label when
#' `evaluate = TRUE`.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param rare_threshold Maximum class proportion used to define rare classes
#' when calculating `rare_recall`.
#' @param scomm_python Optional Python binary used by the `scOMM` backend.
#' If `NULL`, `SCOP_SCOMM_PYTHON` is consulted and reticulate defaults are used
#' otherwise.
#' @param scomm_hidden_nodes,scomm_epochs,scomm_batch_size,scomm_threshold,scomm_seed
#' Parameters passed to the `scOMM` backend.
#'
#' @return A `Seurat` object with `scOMM` predictions stored in metadata and `tools`.
#' @export
#' @examples
#' \dontrun{
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub <- standard_scop(pbmcmultiome_sub, assay = "RNA")
#' ref_cells <- colnames(pbmcmultiome_sub)[1:250]
#' query_cells <- colnames(pbmcmultiome_sub)[251:500]
#' reference <- subset(pbmcmultiome_sub, cells = ref_cells)
#' query <- subset(pbmcmultiome_sub, cells = query_cells)
#' query <- RunscOMM(
#'   srt = query,
#'   reference = reference,
#'   reference_assay = "RNA",
#'   query_assay = "RNA",
#'   reference_label = "CellType",
#'   scomm_epochs = 1
#' )
#' CellDimPlot(
#'   query,
#'   group.by = c(
#'     "CellType",
#'     "scomm_prediction"
#'   ),
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
#'
#' FeatureDimPlot(
#'   query,
#'   features = c(
#'     "scomm_score.B",
#'     "scomm_score.max"
#'   ),
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
#' }
RunscOMM <- function(
  srt,
  reference,
  reference_assay = NULL,
  query_assay = NULL,
  reference_label = NULL,
  features = NULL,
  prediction_prefix = "scomm_",
  evaluate = FALSE,
  truth_col = NULL,
  tool_name = "scOMM",
  rare_threshold = 0.05,
  scomm_python = NULL,
  scomm_hidden_nodes = c(128, 64),
  scomm_epochs = 10,
  scomm_batch_size = 32,
  scomm_threshold = 0.5,
  scomm_seed = 11,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!inherits(reference, "Seurat")) {
    log_message(
      "{.arg reference} is not a {.cls Seurat}",
      message_type = "error"
    )
  }

  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt)
  tool_name <- tool_name %||% "scOMM"
  if (!reference_assay %in% SeuratObject::Assays(reference)) {
    log_message(
      "{.arg reference_assay} not found in {.arg reference}",
      message_type = "error"
    )
  }
  if (!query_assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg query_assay} not found in {.arg srt}",
      message_type = "error"
    )
  }

  scomm_res <- run_scomm(
    reference = reference,
    query = srt,
    reference_assay = reference_assay,
    query_assay = query_assay,
    reference_label = reference_label,
    features = features,
    python = scomm_python,
    hidden_nodes = scomm_hidden_nodes,
    epochs = scomm_epochs,
    batch_size = scomm_batch_size,
    threshold = scomm_threshold,
    seed = scomm_seed,
    verbose = verbose
  )
  pred_res <- add_prediction_meta(
    srt = srt,
    ids = scomm_res$ids,
    probabilities = scomm_res$probabilities,
    prediction_prefix = prediction_prefix
  )
  srt <- pred_res$srt

  if (isTRUE(evaluate)) {
    if (is.null(truth_col) || !truth_col %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg truth_col} must be provided in {.arg srt@meta.data} when {.arg evaluate = TRUE}",
        message_type = "error"
      )
    }
    eval_res <- collect_mapping_metrics(
      srt = srt,
      predicted_col = pred_res$predicted_col,
      truth_col = truth_col,
      probability_col = pred_res$probability_col,
      rare_threshold = rare_threshold
    )
  } else {
    eval_res <- NULL
  }

  srt@tools[[tool_name]] <- list(
    method = "scOMM",
    predicted_col = pred_res$predicted_col,
    probability_col = pred_res$probability_col,
    probability_cols = pred_res$probability_cols,
    truth_col = truth_col,
    metrics = eval_res,
    reference_assay = reference_assay,
    query_assay = query_assay,
    reference_label = reference_label,
    features = scomm_res$features
  )
  srt
}

add_prediction_meta <- function(
  srt,
  ids,
  probabilities = NULL,
  prediction_prefix = "scomm_"
) {
  pred_col <- paste0(prediction_prefix, "prediction")
  srt[[pred_col]] <- factor(ids)
  prob_col <- NULL
  prob_cols <- character(0)
  if (!is.null(probabilities)) {
    prob_df <- as.data.frame(probabilities, check.names = FALSE)
    prob_cols <- paste0(
      prediction_prefix,
      "score.",
      make.names(colnames(prob_df))
    )
    colnames(prob_df) <- prob_cols
    rownames(prob_df) <- colnames(srt)
    srt <- SeuratObject::AddMetaData(srt, metadata = prob_df)
    prob_col <- paste0(prediction_prefix, "score.max")
    check_r("matrixStats", verbose = FALSE)
    srt[[prob_col]] <- matrixStats::rowMaxs(as.matrix(prob_df), na.rm = TRUE)
  }
  list(
    srt = srt,
    predicted_col = pred_col,
    probability_col = prob_col,
    probability_cols = prob_cols
  )
}

run_scomm <- function(
  reference,
  query,
  reference_assay,
  query_assay,
  reference_label,
  features = NULL,
  python = NULL,
  hidden_nodes = c(128, 64),
  epochs = 10,
  batch_size = 32,
  threshold = 0.5,
  seed = 11,
  verbose = TRUE
) {
  python <- python %||% Sys.getenv("SCOP_SCOMM_PYTHON", unset = "")
  if (!nzchar(python)) {
    cached_modules <- tryCatch(
      getOption("scop_env_cache", default = NULL)[["modules"]],
      error = function(...) NULL
    )
    env_modules <- unique(c(cached_modules %||% character(0), "scomm"))
    PrepareEnv(modules = env_modules)
    python <- tryCatch(
      getOption("scop_env_cache", default = NULL)[["python"]],
      error = function(...) NULL
    )
    python <- python %||%
      tryCatch(
        conda_python(envname = get_envname(), conda = resolve_conda("auto")),
        error = function(...) ""
      )
  }
  if (nzchar(python)) {
    configure_python_runtime(python)
  }
  if (
    is.null(reference_label) ||
      !reference_label %in% colnames(reference@meta.data)
  ) {
    log_message(
      "{.arg reference_label} must be present in {.arg reference@meta.data} for {.val method = 'scOMM'}.",
      message_type = "error"
    )
  }
  features <- features %||%
    SeuratObject::VariableFeatures(reference, assay = reference_assay)
  ref_mat <- GetAssayData5(reference, assay = reference_assay, layer = "data")
  query_mat <- GetAssayData5(query, assay = query_assay, layer = "data")
  features <- Reduce(
    intersect,
    list(features, rownames(ref_mat), rownames(query_mat))
  )
  if (length(features) < 10) {
    log_message(
      "Need at least 10 shared features for {.val method = 'scOMM'}",
      message_type = "error"
    )
  }
  ref_mat <- ref_mat[features, , drop = FALSE]
  query_mat <- query_mat[features, , drop = FALSE]
  labels <- factor(reference[[reference_label, drop = TRUE]])
  names(labels) <- colnames(reference)
  split_data <- prepare_scomm_data(
    scale_data = ref_mat,
    labels = labels,
    seed = seed,
    verbose = verbose
  )
  if (!identical(threshold, "AUTO") && nzchar(python)) {
    preds <- run_scomm_subprocess(
      split_data = split_data,
      query_mat = query_mat,
      python = python,
      hidden_nodes = hidden_nodes,
      epochs = epochs,
      batch_size = batch_size,
      threshold = threshold,
      seed = seed
    )
    return(list(
      ids = preds$ids,
      probabilities = preds$probabilities,
      features = features
    ))
  }
  if (is.null(python) || !nzchar(python)) {
    log_message(
      "A Python environment with {.pkg tensorflow} is required for {.val method = 'scOMM'}.",
      message_type = "error"
    )
  }
  tensorflow_ok <- tryCatch(
    {
      reticulate::use_python(python, required = TRUE)
      reticulate::import("tensorflow", delay_load = FALSE)
      TRUE
    },
    error = function(...) FALSE
  )
  if (!isTRUE(tensorflow_ok)) {
    log_message(
      c(
        "Python {.pkg tensorflow} is required for {.val method = 'scOMM'}.",
        "Run {.code PrepareEnv(modules = 'scomm')} or provide {.arg scomm_python} pointing to a Python with tensorflow installed."
      ),
      message_type = "error"
    )
  }
  patch_keras_categorical()
  if (
    requireNamespace("keras", quietly = TRUE) &&
      "py_require_legacy_keras" %in% getNamespaceExports("keras")
  ) {
    tryCatch(
      keras::py_require_legacy_keras(),
      error = function(...) NULL
    )
  }
  reticulate::use_python(python, required = TRUE)
  tryCatch(
    reticulate::py_run_string(
      "import sys\nsys.modules.setdefault('jax', None)\nsys.modules.setdefault('jaxlib', None)\n"
    ),
    error = function(...) NULL
  )
  model <- train_scomm(
    out = split_data,
    hidden_nodes = hidden_nodes,
    epochs = epochs,
    seed = seed,
    batch_size = batch_size,
    verbose = verbose
  )
  preds <- predict_scomm(
    dnn_model = model,
    model.data = split_data,
    query.data = query_mat,
    threshold = threshold
  )
  list(
    ids = preds$ids,
    probabilities = preds$probabilities,
    features = features
  )
}

prepare_scomm_data <- function(
  scale_data,
  labels,
  seed = 11,
  verbose = TRUE
) {
  keep <- !is.na(labels)
  if (!any(keep)) {
    log_message(
      "No valid labels available for {.val method = 'scOMM'} training.",
      message_type = "error"
    )
  }
  labels <- droplevels(factor(labels[keep]))
  train_x <- t(as.matrix(scale_data[, keep, drop = FALSE]))
  classes <- levels(labels)
  if (length(classes) < 2) {
    log_message(
      "{.val method = 'scOMM'} requires at least two reference classes.",
      message_type = "error"
    )
  }
  train_y <- matrix(
    0,
    nrow = nrow(train_x),
    ncol = length(classes),
    dimnames = list(NULL, classes)
  )
  train_y[cbind(seq_len(nrow(train_y)), as.integer(labels))] <- 1
  set.seed(seed)
  log_message("Preparing scOMM training matrices...", verbose = verbose)
  list(
    train_x = train_x,
    train_y = train_y,
    classes = classes
  )
}

run_scomm_subprocess <- function(
  split_data,
  query_mat,
  python,
  hidden_nodes,
  epochs,
  batch_size,
  threshold,
  seed
) {
  tmpdir <- tempfile(pattern = "scop-scomm-")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  if (identical(threshold, "AUTO")) {
    log_message(
      "{.arg threshold = 'AUTO'} is not supported in the external {.pkg python} scOMM runner.",
      message_type = "error"
    )
  }

  train_x <- split_data$train_x
  class_names <- as.character(split_data$classes)
  query_x <- t(as.matrix(query_mat))
  missing_features <- setdiff(colnames(train_x), colnames(query_x))
  if (length(missing_features) > 0) {
    padding <- matrix(
      0,
      nrow = nrow(query_x),
      ncol = length(missing_features),
      dimnames = list(NULL, missing_features)
    )
    query_x <- cbind(query_x, padding)
  }
  query_x <- as.matrix(query_x[, colnames(train_x), drop = FALSE])
  threshold_vec <- if (length(threshold) == 1L) {
    rep(as.numeric(threshold), length(class_names))
  } else {
    threshold[class_names]
  }
  names(threshold_vec) <- class_names

  train_x_file <- file.path(tmpdir, "train_x.csv")
  train_y_file <- file.path(tmpdir, "train_y.csv")
  query_x_file <- file.path(tmpdir, "query_x.csv")
  prob_file <- file.path(tmpdir, "probabilities.csv")
  class_file <- file.path(tmpdir, "classes.txt")
  threshold_file <- file.path(tmpdir, "thresholds.txt")
  script_file <- file.path(tmpdir, "run_scomm_child.py")

  utils::write.table(
    train_x,
    train_x_file,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  utils::write.table(
    as.matrix(split_data$train_y),
    train_y_file,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  utils::write.table(
    query_x,
    query_x_file,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  writeLines(class_names, class_file, useBytes = TRUE)
  writeLines(as.character(threshold_vec), threshold_file, useBytes = TRUE)
  writeLines(scomm_python_script(), con = script_file, useBytes = TRUE)

  args <- c(
    "-i",
    scomm_subprocess_env(python),
    python,
    script_file,
    train_x_file,
    train_y_file,
    query_x_file,
    prob_file,
    as.character(length(class_names)),
    paste(hidden_nodes, collapse = ","),
    as.character(epochs),
    as.character(batch_size),
    as.character(seed)
  )
  output <- suppressWarnings(
    system2(
      command = "/usr/bin/env",
      args = shQuote(args),
      stdout = TRUE,
      stderr = TRUE
    )
  )
  status <- attr(output, "status") %||% 0L
  if (!identical(as.integer(status), 0L) || !file.exists(prob_file)) {
    log_message(
      c(
        "{.val method = 'scOMM'} subprocess failed.",
        output
      ),
      message_type = "error"
    )
  }
  probs <- as.matrix(utils::read.csv(
    prob_file,
    header = FALSE,
    check.names = FALSE
  ))
  colnames(probs) <- class_names
  predicted <- scomm_assign_labels(probs, threshold_vec)
  list(
    ids = factor(predicted),
    probabilities = as.data.frame(
      probs,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

scomm_python_script <- function() {
  c(
    "import os",
    "import sys",
    "import numpy as np",
    "os.environ.setdefault('OMP_NUM_THREADS', '1')",
    "os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')",
    "os.environ.setdefault('MKL_NUM_THREADS', '1')",
    "os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')",
    "os.environ.setdefault('NUMEXPR_NUM_THREADS', '1')",
    "os.environ.setdefault('NUMBA_NUM_THREADS', '1')",
    "import tensorflow as tf",
    "train_x = np.loadtxt(sys.argv[1], delimiter=',', ndmin=2).astype('float32')",
    "train_y = np.loadtxt(sys.argv[2], delimiter=',', ndmin=2).astype('float32')",
    "query_x = np.loadtxt(sys.argv[3], delimiter=',', ndmin=2).astype('float32')",
    "prob_file = sys.argv[4]",
    "n_classes = int(sys.argv[5])",
    "hidden_nodes = [int(x) for x in sys.argv[6].split(',') if x]",
    "epochs = int(sys.argv[7])",
    "batch_size = int(sys.argv[8])",
    "seed = int(sys.argv[9])",
    "tf.keras.utils.set_random_seed(seed)",
    "model = tf.keras.Sequential()",
    "model.add(tf.keras.layers.Input(shape=(train_x.shape[1],)))",
    "for units in hidden_nodes:",
    "    model.add(tf.keras.layers.Dense(units=units, activation='relu'))",
    "    model.add(tf.keras.layers.Dropout(rate=0.2))",
    "model.add(tf.keras.layers.Dense(units=n_classes, activation='softmax'))",
    "model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), loss='categorical_crossentropy', metrics=['accuracy'])",
    "model.fit(train_x, train_y, epochs=epochs, batch_size=batch_size, validation_split=0.2, verbose=0)",
    "probs = model.predict(query_x, verbose=0)",
    "np.savetxt(prob_file, probs, delimiter=',')"
  )
}

scomm_assign_labels <- function(
  probabilities,
  threshold_vec,
  unclassified = "unclassified"
) {
  probabilities <- as.matrix(probabilities)
  thresholds <- threshold_vec[colnames(probabilities)]
  if (anyNA(probabilities) || anyNA(thresholds)) {
    return(apply(probabilities, 1, function(x) {
      valid <- x >= thresholds[names(x)]
      if (any(valid)) {
        names(which.max(x[valid]))[[1]]
      } else {
        unclassified
      }
    }))
  }
  valid <- sweep(probabilities, 2, thresholds, FUN = ">=")
  has_valid <- rowSums(valid) > 0L
  out <- rep(unclassified, nrow(probabilities))
  if (any(has_valid)) {
    masked <- probabilities[has_valid, , drop = FALSE]
    masked[!valid[has_valid, , drop = FALSE]] <- -Inf
    out[has_valid] <- colnames(probabilities)[
      max.col(masked, ties.method = "first")
    ]
  }
  out
}

scomm_subprocess_env <- function(python) {
  python <- normalizePath(python, mustWork = FALSE)
  python_dir <- dirname(python)
  env_path <- dirname(python_dir)
  join_env <- function(var, values) {
    values <- values[nzchar(values)]
    values <- values[file.exists(values)]
    current <- Sys.getenv(var, unset = "")
    current_values <- strsplit(current, .Platform$path.sep, fixed = TRUE)[[1]]
    current_values <- current_values[nzchar(current_values)]
    paste(unique(c(values, current_values)), collapse = .Platform$path.sep)
  }
  env <- c(
    HOME = Sys.getenv("HOME", unset = "~"),
    TMPDIR = Sys.getenv("TMPDIR", unset = tempdir()),
    LANG = Sys.getenv("LANG", unset = "C.UTF-8"),
    LC_ALL = Sys.getenv("LC_ALL", unset = "C.UTF-8"),
    RETICULATE_PYTHON = python,
    RETICULATE_PYTHON_ENV = "",
    PYTHONNOUSERSITE = "1",
    PYTHONPATH = "",
    PYTHONHOME = "",
    PIP_USER = "0",
    PATH = join_env(
      "PATH",
      c(
        python_dir,
        file.path(env_path, "bin"),
        file.path(env_path, "Library", "bin"),
        file.path(env_path, "Scripts")
      )
    ),
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    NUMBA_NUM_THREADS = "1",
    KMP_WARNINGS = "0",
    KMP_DUPLICATE_LIB_OK = "TRUE"
  )
  if (is_linux()) {
    env <- c(
      env,
      LD_LIBRARY_PATH = join_env(
        "LD_LIBRARY_PATH",
        c(
          file.path(env_path, "lib"),
          file.path(env_path, "lib64"),
          file.path(env_path, "Library", "lib")
        )
      )
    )
  }
  if (is_osx()) {
    env <- c(
      env,
      DYLD_FALLBACK_LIBRARY_PATH = join_env(
        "DYLD_FALLBACK_LIBRARY_PATH",
        c(
          file.path(env_path, "lib"),
          file.path(env_path, "lib64")
        )
      )
    )
  }
  unname(paste(names(env), env, sep = "="))
}

patch_keras_categorical <- function() {
  if (!requireNamespace("keras", quietly = TRUE)) {
    return(invisible(NULL))
  }
  ns <- asNamespace("keras")
  to_categorical_compat <- function(
    y = NULL,
    x = y,
    num_classes = NULL,
    dtype = "float32"
  ) {
    vec <- as.integer(y %||% x)
    keep <- !is.na(vec)
    vec_use <- vec[keep]
    if (length(vec_use) == 0) {
      return(matrix(numeric(0), nrow = 0))
    }
    if (is.null(num_classes)) {
      num_classes <- max(vec_use)
    }
    tf <- reticulate::import("tensorflow", delay_load = FALSE)
    res <- tryCatch(
      tf$keras$utils$to_categorical(
        as.integer(vec_use - 1L),
        num_classes = as.integer(num_classes),
        dtype = dtype
      ),
      error = function(...) {
        tf$keras$utils$to_categorical(
          as.integer(vec_use - 1L),
          num_classes = as.integer(num_classes)
        )
      }
    )
    out <- reticulate::py_to_r(res)
    if (!all(keep)) {
      full_out <- matrix(0, nrow = length(vec), ncol = ncol(out))
      full_out[keep, ] <- out
      return(full_out)
    }
    out
  }
  env_candidates <- unique(list(
    ns,
    parent.env(ns)
  ))
  for (env_i in env_candidates) {
    if (!exists("to_categorical", envir = env_i, inherits = FALSE)) {
      next
    }
    tryCatch(
      {
        get("unlockBinding", envir = baseenv())("to_categorical", env_i)
        assign("to_categorical", to_categorical_compat, envir = env_i)
        get("lockBinding", envir = baseenv())("to_categorical", env_i)
      },
      error = function(...) NULL
    )
  }
  invisible(NULL)
}

train_scomm <- function(
  out,
  hidden_nodes,
  epochs = 10,
  seed = 11,
  batch_size = 32,
  activation = "relu",
  add_dropout = TRUE,
  pct_dropout = 0.2,
  lr = 0.001,
  l1 = 0,
  l2 = 0,
  verbose = TRUE
) {
  tf <- reticulate::import("tensorflow", delay_load = FALSE)
  train_x <- as.matrix(out$train_x)
  train_y <- as.matrix(out$train_y)
  input_dim <- ncol(train_x)
  output_dim <- ncol(train_y)

  tryCatch(
    tf$keras$utils$set_random_seed(as.integer(seed)),
    error = function(...) NULL
  )

  model <- tf$keras$Sequential()
  model$add(
    tf$keras$layers$Input(
      shape = reticulate::tuple(as.integer(input_dim))
    )
  )

  reg <- if (isTRUE(l1 > 0 || l2 > 0)) {
    tf$keras$regularizers$l1_l2(
      l1 = as.numeric(l1),
      l2 = as.numeric(l2)
    )
  } else {
    NULL
  }

  for (units in hidden_nodes) {
    model$add(
      tf$keras$layers$Dense(
        units = as.integer(units),
        activation = activation,
        kernel_regularizer = reg
      )
    )
    if (isTRUE(add_dropout) && pct_dropout > 0) {
      model$add(
        tf$keras$layers$Dropout(rate = as.numeric(pct_dropout))
      )
    }
  }

  model$add(
    tf$keras$layers$Dense(
      units = as.integer(output_dim),
      activation = "softmax"
    )
  )

  model$compile(
    optimizer = tf$keras$optimizers$Adam(learning_rate = as.numeric(lr)),
    loss = "categorical_crossentropy",
    metrics = list("accuracy")
  )

  model$fit(
    x = train_x,
    y = train_y,
    epochs = as.integer(epochs),
    batch_size = as.integer(batch_size),
    validation_split = 0.2,
    verbose = if (isTRUE(verbose)) 2L else 0L
  )

  model
}

predict_scomm <- function(
  dnn_model,
  model.data,
  query.data,
  threshold = 0.5
) {
  train_x <- model.data$train_x
  features <- colnames(train_x)
  test_x <- t(as.matrix(query.data))
  do.genes <- setdiff(features, colnames(test_x))
  if (length(do.genes) > 0) {
    do.genes_x <- matrix(0, ncol = length(do.genes), nrow = nrow(test_x))
    colnames(do.genes_x) <- do.genes
    test_x <- cbind(test_x, do.genes_x)
    test_x <- test_x[, colnames(train_x), drop = FALSE]
  } else {
    test_x <- as.matrix(test_x[, features, drop = FALSE])
  }

  probs <- tryCatch(
    stats::predict(dnn_model, test_x, verbose = 0),
    error = function(...) dnn_model$predict(test_x, verbose = as.integer(0))
  )
  probs <- tryCatch(reticulate::py_to_r(probs), error = function(...) probs)
  probs <- as.data.frame(probs, check.names = FALSE, stringsAsFactors = FALSE)
  class_names <- as.character(model.data$classes)
  if (ncol(probs) == length(class_names) + 1L) {
    colnames(probs) <- c("unclassified", class_names)
    probs_use <- probs[, class_names, drop = FALSE]
  } else if (ncol(probs) == length(class_names)) {
    colnames(probs) <- class_names
    probs_use <- probs
  } else {
    log_message(
      "{.val method = 'scOMM'} returned an unexpected number of probability columns.",
      message_type = "error"
    )
  }

  if (identical(threshold, "AUTO")) {
    check_r("mereulab/scOMM", verbose = FALSE)
    threshold_vec <- get_namespace_fun("scOMM", "find_optimal_thresholds")(
      dnn_model = dnn_model,
      model.data = model.data
    )
  } else if (length(threshold) == 1L) {
    threshold_vec <- rep(as.numeric(threshold), length(class_names))
    names(threshold_vec) <- class_names
  } else {
    threshold_vec <- threshold
    if (!all(class_names %in% names(threshold_vec))) {
      log_message(
        "Named {.arg threshold} must include all scOMM classes.",
        message_type = "error"
      )
    }
    threshold_vec <- threshold_vec[class_names]
  }

  predicted <- scomm_assign_labels(probs_use, threshold_vec)

  list(
    ids = factor(predicted),
    probabilities = probs_use
  )
}
