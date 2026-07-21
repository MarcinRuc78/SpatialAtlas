#!/usr/bin/env Rscript

options(warn = 1)
suppressPackageStartupMessages({
  library(hdf5r)
  library(Matrix)
})

input_h5ad <- Sys.getenv(
  "ADRENAL_HD_H5AD",
  "/package/examples/02_visium_hd_8um/output/sample_adrenal_female_hd_8um.h5ad"
)
output_rds <- Sys.getenv(
  "ADRENAL_HD_SHINY_RDS",
  "/work/source_data/adrenal_visium_hd_shiny.rds"
)

read_vector <- function(file, path) {
  file[[path]]$read()
}

file <- H5File$new(input_h5ad, mode = "r")
on.exit(file$close_all(), add = TRUE)

data <- read_vector(file, "X/data")
indices <- as.integer(read_vector(file, "X/indices"))
indptr <- as.integer(read_vector(file, "X/indptr"))
genes <- as.character(read_vector(file, "var/_index"))
barcodes <- as.character(read_vector(file, "obs/_index"))
coordinates <- file[["obsm/spatial"]]$read()

# AnnData stores X as observation-by-gene CSR. The same arrays are a
# gene-by-observation dgCMatrix after swapping the dimensions, without a dense
# intermediate or a numerical transformation.
expression <- new(
  "dgCMatrix",
  i = indices,
  p = indptr,
  x = as.numeric(data),
  Dim = as.integer(c(length(genes), length(barcodes))),
  Dimnames = list(genes, barcodes)
)

if (!identical(dim(coordinates), c(length(barcodes), 2L))) {
  coordinates <- t(coordinates)
}
stopifnot(
  identical(dim(expression), c(length(genes), length(barcodes))),
  identical(dim(coordinates), c(length(barcodes), 2L)),
  all(is.finite(coordinates)),
  "Cyp11b2" %in% genes
)

object <- list(
  expression = expression,
  coordinates = unname(coordinates),
  genes = genes,
  barcodes = barcodes,
  coordinate_scale = 0.205592,
  source_shape = c(profiles = length(barcodes), genes = length(genes)),
  source_format = "AnnData CSR; log-normalized expression retained without transformation"
)
class(object) <- c("adrenal_visium_hd_shiny", "list")

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(object, output_rds, compress = "xz")
cat(
  sprintf(
    "Saved %s: %s profiles x %s genes; %.1f MiB\n",
    output_rds,
    format(length(barcodes), big.mark = ","),
    format(length(genes), big.mark = ","),
    file.info(output_rds)$size / 2^20
  )
)
