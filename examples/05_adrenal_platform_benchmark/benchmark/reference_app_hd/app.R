suppressPackageStartupMessages({
  library(Cairo)
  library(ggplot2)
  library(grid)
  library(Matrix)
  library(png)
  library(shiny)
})

object <- readRDS(Sys.getenv(
  "SHINY_HD_RDS",
  "/work/source_data/adrenal_visium_hd_shiny.rds"
))
histology <- png::readPNG(Sys.getenv(
  "SHINY_HD_IMAGE",
  "/data/tissue_adrenal_female_hd_8um.png"
), native = TRUE)

stopifnot(
  inherits(object$expression, "dgCMatrix"),
  ncol(object$expression) == nrow(object$coordinates),
  "Cyp11b2" %in% object$genes
)

features <- sort(object$genes)
gene_index <- setNames(seq_along(object$genes), object$genes)
image_width <- ncol(histology)
image_height <- nrow(histology)
coordinates <- data.frame(
  x = object$coordinates[, 1] * object$coordinate_scale,
  y = image_height - object$coordinates[, 2] * object$coordinate_scale
)

get_expression <- function(gene) {
  as.numeric(object$expression[gene_index[[gene]], , drop = TRUE])
}

build_spatial_plot <- function(gene, point_size = 0.34, minimum_alpha = 0.10,
                               expression = NULL) {
  if (is.null(expression)) {
    expression <- get_expression(gene)
  }
  plot_data <- coordinates
  plot_data$expression <- expression
  ggplot(plot_data, aes(x = x, y = y)) +
    annotation_raster(
      histology,
      xmin = 0,
      xmax = image_width,
      ymin = 0,
      ymax = image_height
    ) +
    geom_point(
      aes(fill = expression, alpha = expression),
      shape = 21,
      colour = "transparent",
      size = point_size,
      stroke = 0
    ) +
    scale_fill_gradientn(
      colours = grDevices::hcl.colors(100, palette = "Inferno"),
      name = gene
    ) +
    scale_alpha(range = c(minimum_alpha, 0.92), guide = "none") +
    coord_fixed(
      xlim = c(0, image_width),
      ylim = c(0, image_height),
      expand = FALSE
    ) +
    theme_void(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

ui <- fluidPage(
  titlePanel("Mouse adrenal gland — Visium HD"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "feature",
        "Feature to visualize",
        choices = features,
        selected = "Cyp11b2",
        selectize = FALSE,
        width = "100%"
      ),
      sliderInput(
        "point_size",
        "Point size",
        min = 0.1,
        max = 1.2,
        value = 0.34,
        step = 0.02,
        width = "100%"
      ),
      sliderInput(
        "minimum_alpha",
        "Minimum alpha",
        min = 0,
        max = 0.5,
        value = 0.10,
        step = 0.01,
        width = "100%"
      )
    ),
    mainPanel(plotOutput("spatial_plot", width = "100%", height = "820px"))
  )
)

server <- function(input, output, session) {
  output$spatial_plot <- renderPlot(
    build_spatial_plot(input$feature, input$point_size, input$minimum_alpha),
    width = 900,
    height = 900,
    res = 96,
    bg = "white"
  )
}

shinyApp(ui = ui, server = server)
