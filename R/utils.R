kegg_get <- function(url) {
  temp <- tempfile()
  on.exit(unlink(temp))
  download(url = url, destfile = temp)
  content <- as.data.frame(
    do.call(
      rbind,
      strsplit(readLines(temp), split = "\t")
    )
  )
  content
}

col2hex <- function(cname) {
  col_mat <- grDevices::col2rgb(cname)
  grDevices::rgb(
    red = col_mat[1, ] / 255,
    green = col_mat[2, ] / 255,
    blue = col_mat[3, ] / 255
  )
}

select_cells <- function(obj, celltypes, group.by) {
  metadata <- obj@meta.data
  cells_c <- c()
  for (celltype in celltypes) {
    cells_c <- c(
      cells_c,
      rownames(metadata[metadata[[group.by]] == celltype, ])
    )
  }
  cells_c
}
