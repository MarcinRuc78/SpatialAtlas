#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Cairo)
  library(Matrix)
  library(Seurat)
  library(grid)
})

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[1]), ".."), mustWork = FALSE)
args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "/work/source_data/adrenal_visium_seurat.rds"
output_dir <- if (length(args) >= 2) args[[2]] else "/work/prepared/intermediate"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

adrenal <- readRDS(input_rds)
assay_name <- DefaultAssay(adrenal)
expression <- GetAssayData(adrenal, assay = assay_name, layer = "data")
cells <- colnames(expression)

Matrix::writeMM(t(expression), file.path(output_dir, "expression_log1p.mtx"))
write.table(
  rownames(expression), file.path(output_dir, "features.tsv"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
  cells, file.path(output_dir, "barcodes.tsv"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

coordinates <- GetTissueCoordinates(adrenal@images$slice1, scale = "lowres")
coordinates <- coordinates[cells, c("imagecol", "imagerow"), drop = FALSE]
coordinates$barcode <- rownames(coordinates)
write.csv(coordinates[, c("barcode", "imagecol", "imagerow")],
          file.path(output_dir, "coordinates_lowres.csv"), row.names = FALSE)

image_grob <- GetImage(adrenal@images$slice1)
image_width <- ncol(image_grob$raster)
image_height <- nrow(image_grob$raster)
CairoPNG(file.path(output_dir, "tissue_adrenal_visium.png"),
         width = image_width, height = image_height, bg = "white")
grid.newpage()
grid.draw(image_grob)
dev.off()

identities <- data.frame(
  barcode = cells,
  cluster = as.character(Idents(adrenal)[cells]),
  stringsAsFactors = FALSE
)
write.csv(identities, file.path(output_dir, "identities.csv"), row.names = FALSE)

cat(sprintf("profiles=%d\n", ncol(expression)))
cat(sprintf("genes=%d\n", nrow(expression)))
cat(sprintf("image_width=%d\n", image_width))
cat(sprintf("image_height=%d\n", image_height))
cat(sprintf("assay=%s\n", assay_name))
