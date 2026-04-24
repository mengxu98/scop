metric_entropy <- function(x) {
  x <- x[x > 0]
  if (length(x) == 0) {
    return(0)
  }
  p <- x / sum(x)
  -sum(p * log(p))
}
