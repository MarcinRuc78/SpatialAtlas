#!/usr/bin/env Rscript

options(warn = 1)
app_env <- new.env(parent = globalenv())
source("/work/benchmark/reference_app_hd/app.R", local = app_env)

genes <- c("Npy", "Cyp11b1", "Cyp11b2", "Th")
repetitions <- 8L
rows <- list()

measure_gene <- function(gene, iteration, keep_row = TRUE) {
  t0 <- proc.time()[["elapsed"]]
  expression <- app_env$get_expression(gene)
  extraction_s <- proc.time()[["elapsed"]] - t0

  t1 <- proc.time()[["elapsed"]]
  plot <- app_env$build_spatial_plot(
    gene,
    point_size = 0.34,
    minimum_alpha = 0.10,
    expression = expression
  )
  figure_build_s <- proc.time()[["elapsed"]] - t1

  output_file <- tempfile(fileext = ".png")
  t2 <- proc.time()[["elapsed"]]
  CairoPNG(filename = output_file, width = 900, height = 900, res = 96)
  print(plot)
  dev.off()
  response_render_s <- proc.time()[["elapsed"]] - t2
  response_bytes <- file.info(output_file)$size
  unlink(output_file)

  if (keep_row) {
    data.frame(
      implementation = "seurat_shiny_visium_hd",
      gene = gene,
      iteration = iteration,
      extraction_s = extraction_s,
      figure_build_s = figure_build_s,
      response_render_s = response_render_s,
      total_server_s = extraction_s + figure_build_s + response_render_s,
      response_bytes = response_bytes
    )
  }
}

for (gene in genes) measure_gene(gene, 0L, keep_row = FALSE)
for (iteration in seq_len(repetitions)) {
  for (gene in genes) {
    rows[[length(rows) + 1L]] <- measure_gene(gene, iteration)
  }
}

result <- do.call(rbind, rows)
write.csv(result, "/work/benchmark/shiny_hd_core_results.csv", row.names = FALSE)
print(aggregate(
  result[c("extraction_s", "figure_build_s", "response_render_s", "total_server_s", "response_bytes")],
  by = list(implementation = result$implementation),
  FUN = median
))
