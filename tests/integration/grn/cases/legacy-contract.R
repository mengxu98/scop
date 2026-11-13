run_case <- function(output) {
  set.seed(1201)
  expression <- matrix(stats::rpois(12L * 48L, 3), 12L,
    dimnames = list(paste0("G", seq_len(12L)), paste0("C", seq_len(48L))))
  regulators <- rownames(expression)[1:4]
  targets <- rownames(expression)[5:10]
  direct <- scop::RunGRNBoost2(expression, regulators = regulators, targets = targets,
    genes_in = "rows", n_rounds = 20, verbose = FALSE)
  edges <- scop::RunGRN(expression, regulators = regulators, targets = targets,
    methods = "grnboost2", n_rounds = 20, verbose = FALSE)
  stopifnot(identical(edges, direct))
  result <- scop::RunGRN(expression, regulators = regulators, targets = targets,
    methods = "grnboost2", method_args = list(n_rounds = 20),
    return_type = "result", work_dir = file.path(output, "work"), verbose = FALSE)
  stopifnot(identical(scop::GetGRN(result, what = "native"), direct))
  result
}
