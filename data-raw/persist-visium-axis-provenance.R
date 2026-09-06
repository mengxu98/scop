# Run from the package root. Migrate the verified, existing image attribute
# into the Seurat object misc slot without changing measurements or pixels.
path <- "data/visium_human_pancreas_sub.rda"
load(path)
before <- visium_human_pancreas_sub
stopifnot(identical(attr(before@images$slice1, "coords_x_orientation"), "horizontal"))
visium_human_pancreas_sub@misc$spatial_image_axes$slice1 <- "horizontal"
restored <- visium_human_pancreas_sub
restored@misc <- before@misc
stopifnot(identical(restored, before))
save(visium_human_pancreas_sub, file = path, compress = "xz")
