#!/usr/bin/env Rscript

options(warn = 1)
suppressPackageStartupMessages({
  library(hdf5r)
  library(jsonlite)
  library(Matrix)
})

input_h5ad <- Sys.getenv(
  "ADRENAL_HD_H5AD",
  "/package/examples/02_visium_hd_8um/output/sample_adrenal_female_hd_8um.h5ad"
)
input_rds <- Sys.getenv(
  "ADRENAL_HD_SHINY_RDS",
  "/work/source_data/adrenal_visium_hd_shiny.rds"
)
output_json <- Sys.getenv(
  "ADRENAL_HD_EQUIVALENCE_JSON",
  "/work/benchmark/hd_object_equivalence.json"
)

file <- H5File$new(input_h5ad, mode = "r")
on.exit(file$close_all(), add = TRUE)
object <- readRDS(input_rds)

h5_data <- as.numeric(file[["X/data"]]$read())
h5_indices <- as.integer(file[["X/indices"]]$read())
h5_indptr <- as.integer(file[["X/indptr"]]$read())
h5_genes <- as.character(file[["var/_index"]]$read())
h5_barcodes <- as.character(file[["obs/_index"]]$read())
h5_coordinates <- file[["obsm/spatial"]]$read()
if (!identical(dim(h5_coordinates), c(length(h5_barcodes), 2L))) {
  h5_coordinates <- t(h5_coordinates)
}

checks <- list(
  dimensions_identical = identical(
    dim(object$expression),
    c(length(h5_genes), length(h5_barcodes))
  ),
  sparse_values_identical = identical(object$expression@x, h5_data),
  sparse_indices_identical = identical(object$expression@i, h5_indices),
  sparse_pointers_identical = identical(object$expression@p, h5_indptr),
  genes_identical = identical(object$genes, h5_genes),
  barcodes_identical = identical(object$barcodes, h5_barcodes),
  coordinates_identical = identical(object$coordinates, unname(h5_coordinates))
)
stopifnot(all(unlist(checks, use.names = FALSE)))

markers <- c("Npy", "Cyp11b1", "Cyp11b2", "Th")
marker_summary <- lapply(markers, function(gene) {
  values <- as.numeric(object$expression[match(gene, object$genes), ])
  list(
    gene = gene,
    nonzero_profiles = unname(sum(values != 0)),
    maximum = unname(max(values))
  )
})

payload <- list(
  validated = TRUE,
  profiles = length(h5_barcodes),
  genes = length(h5_genes),
  nonzero_matrix_entries = length(h5_data),
  checks = checks,
  markers = marker_summary
)
dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
write_json(payload, output_json, pretty = TRUE, auto_unbox = TRUE, digits = 15)
cat("HD H5AD-to-RDS equivalence: PASS\n")
cat(sprintf("Wrote %s\n", output_json))
