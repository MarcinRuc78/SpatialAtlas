library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(shiny)
library(miniUI)
library(grid)
library(Cairo)

# Load the Seurat object
object <- readRDS(Sys.getenv(
  "SHINY_SEURAT_RDS",
  "/work/source_data/adrenal_visium_seurat.rds"
))
feature <- "Npy"
slot <- 'data'
alpha <- c(0.1, 1.0)

# Precompute color palettes
SpatialColors <- colorRampPalette(colors = rev(brewer.pal(n = 11, name = "Spectral")))

FeaturePalettes <- list(
  'Spatial' = SpatialColors(n = 100),
  'Seurat' = c('white', 'red')
)

# Get initial assay, features, and coordinates
assay.keys <- Key(object)[Seurat::Assays(object)]
keyed <- sapply(assay.keys, grepl, x = feature)
assay <- if (any(keyed)) names(which(keyed))[1] else DefaultAssay(object)
features <- sort(rownames(GetAssayData(object, layer = slot, assay = assay)))
coords <- GetTissueCoordinates(object@images$slice1)
cells.use <- Cells(object)

# Fetch initial feature data for plot
feature.data <- FetchData(object, vars = feature, cells = cells.use, layer = slot)
plot.data <- cbind(coords, feature.data)

# Define UI
ui <- miniPage(
  titlePanel("Spatial Transcriptomics Atlas of the Mouse Adrenal Gland"),
  miniContentPanel(
    fillRow(
      sidebarPanel(
        sliderInput('alpha', 'Alpha intensity', min = 0, max = max(alpha), value = min(alpha), step = 0.01, width = '100%'),
        sliderInput('pt.size', 'Point size', min = 0, max = 5, value = 1.1, step = 0.1, width = '100%'),
        selectInput('feature', 'Feature to visualize', choices = features, selected = feature, selectize = FALSE, width = '100%'),
        selectInput('palette', 'Color scheme', choices = names(FeaturePalettes), selected = 'Spatial', selectize = FALSE, width = '100%'),
        width = '100%'
      ),
      plotOutput('plot', height = '100%'),
      flex = c(1, 4)
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  plot.env <- reactiveValues(data = plot.data, feature = feature, palette = 'Spatial')
  
  # Observe changes in feature input
  observeEvent(input$feature, {
    feature.use <- input$feature
    feature.data <- FetchData(object, vars = feature.use, cells = cells.use, layer = slot)
    colnames(feature.data) <- feature.use
    plot.env$data <- cbind(coords, feature.data)
    plot.env$feature <- feature.use
  })
  
  observeEvent(input$palette, {
    plot.env$palette <- input$palette
  })
  
  output$plot <- renderPlot({
    SingleSpatialPlot(
      data = plot.env$data,
      image = object@images$slice1,
      col.by = plot.env$feature,
      pt.size.factor = input$pt.size,
      crop = TRUE,
      alpha.by = plot.env$feature
    ) +
      scale_fill_gradientn(name = plot.env$feature, colours = FeaturePalettes[[plot.env$palette]]) +
      theme(legend.position = 'top') +
      scale_alpha(range = c(input$alpha, 1)) +
      guides(alpha = FALSE)
  })
}

# Define SingleSpatialPlot function
SingleSpatialPlot <- function(
    data,
    image,
    cols = NULL,
    image.alpha = 1,
    pt.alpha = NULL,
    crop = TRUE,
    pt.size.factor = NULL,
    stroke = 0.25,
    col.by = NULL,
    alpha.by = NULL,
    cells.highlight = NULL,
    cols.highlight = c('#DE2D26', 'grey50'),
    geom = 'spatial',
    na.value = 'grey50'
) {
  if (!is.null(col.by) && !col.by %in% colnames(data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  col.by <- col.by %iff% paste0("`", col.by, "`")
  alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
  
  plot <- ggplot(data, aes_string(x = colnames(data)[2], y = colnames(data)[1], fill = col.by, alpha = alpha.by)) +
    geom_spatial(
      point.size.factor = pt.size.factor,
      data = data,
      image = image,
      image.alpha = image.alpha,
      crop = crop,
      stroke = stroke
    ) +
    coord_fixed() +
    theme(aspect.ratio = 1)
  
  if (!is.null(cols) && is.null(cells.highlight)) {
    scale <- if (length(cols) == 1 && cols %in% rownames(brewer.pal.info)) {
      scale_fill_brewer(palette = cols, na.value = na.value)
    } else {
      scale_fill_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  
  plot + NoAxes() + theme(panel.background = element_blank())
}

# Define the GeomSpatial ggproto object
GeomSpatial <- ggproto(
  "GeomSpatial",
  Geom,
  required_aes = c("x", "y"),
  default_aes = aes(
    shape = 21,
    colour = "black",
    point.size.factor = 1.0,
    fill = NA,
    alpha = NA,
    stroke = 0.25
  ),
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    data$y <- max(data$y) - data$y + min(data$y)
    data
  },
  draw_key = draw_key_point,
  draw_panel = function(data, panel_scales, coord, image, image.alpha, crop) {
    if (!crop) {
      y.transform <- c(0, nrow(image)) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$x$continuous_range <- c(0, ncol(image))
      panel_scales$y$continuous_range <- c(0, nrow(image))
      panel_scales$y.range <- c(0, nrow(image))
      panel_scales$x.range <- c(0, ncol(image))
    }
    z <- coord$transform(data.frame(x = c(0, ncol(image)), y = c(0, nrow(image))), panel_scales)
    z$y <- -rev(z$y) + 1
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]
    vp <- viewport(
      x = unit(z$x[1], "npc"),
      y = unit(z$y[1], "npc"),
      width = unit(wdth, "npc"),
      height = unit(hgth, "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)
    
    img <- editGrob(grob = img.grob, vp = vp)
    spot.size <- Radius(object = image)
    coords <- coord$transform(data, panel_scales)
    pts <- pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = data$shape,
      size = unit(spot.size, "npc") * data$point.size.factor,
      gp = gpar(
        col = alpha(colour = coords$colour, alpha = coords$alpha),
        fill = alpha(colour = coords$fill, alpha = coords$alpha),
        lwd = coords$stroke)
    )
    vp <- viewport()
    gt <- gTree(vp = vp)
    if (image.alpha > 0) {
      if (image.alpha != 1) {
        img$raster = as.raster(
          x = matrix(
            data = alpha(colour = img$raster, alpha = image.alpha),
            nrow = nrow(img$raster),
            ncol = ncol(img$raster),
            byrow = TRUE)
        )
      }
      gt <- addGrob(gt, child = img)
    }
    gt <- addGrob(gt, child = pts)
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
    return(gt)
  }
)

# Register geom_spatial function
geom_spatial <- function(
    mapping = NULL,
    data = NULL,
    image = NULL,
    image.alpha = 1,
    crop = TRUE,
    stat = "identity",
    position = "identity",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE,
    ...
) {
  layer(
    geom = GeomSpatial,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, image = image, image.alpha = image.alpha, crop = crop, ...)
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
