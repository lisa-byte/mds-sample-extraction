################################################################################
# Shiny App: MDS Sample/Cluster Extraction
# Author: Lisa Weidlich
# Date: 2025-07-07
#
# Description:
# Interactive Shiny app to manually select sample groups on MDS plots computed
# from pairwise genetic distance matrices. Draw polygons; points inside (or on
# the boundary) are selected. Selected Taxon IDs are printed and saved as
# groupN.txt files in the output directory.
#
# Minimal input:
# (1) Distance matrix: TASSEL IBS format
# (2) Key file (tab-delimited): columns Taxon, <group>, <group>_color (hex);
#     optional: 'proportion missing' in [0,1] used with `propmis` to filter.
# (3) Output directory (wd)
#
# Quick start:
# 1) Edit the USER CONFIG block below.
# 2) Run app.
# 3) Draw polygons (choose polygon symbol before drawing) to export groups.
################################################################################

# ============================ Packages ========================================

# List of required packages
packages <- c(
  "data.table", "shiny", "leaflet", "leaflet.extras", "sf",
  "R.utils", "RColorBrewer", "Hmisc"
)

# Install any that are not present
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
}

# Load them all
lapply(packages, library, character.only = TRUE)

# ============================ User Config =====================================

PERSONAL_DIRECTORY <- "/enter_directory/example_data_region"  # <<< EDIT (obligatory)
propmis <- 0.8             # keep samples with 'proportion missing' < propmis
color_group <- "region"        # key column for colors; needs "<group>_color"

# derived paths
wd              <- file.path(PERSONAL_DIRECTORY, "output")
distance_matrix <- file.path(PERSONAL_DIRECTORY, "data", "example_IBS.dist.txt")
MDS_KEY         <- file.path(PERSONAL_DIRECTORY, "data", "example_key.tsv")

# ============================ MDS (TRG minimal) ===============================
# Credit: K. Swarts (2019); used with permission. Minimal subset adapted for app.

uniqueCols <- function(N, opacity = "80") {
  if (N < 75) {
    qual <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
    col_max <- unlist(mapply(RColorBrewer::brewer.pal, qual$maxcolors, rownames(qual)))
    paste(sample(col_max, N), opacity, sep = "")
  } else {
    paste(substring(text = rainbow(n = N, alpha = 1), 1, 7), opacity, sep = "")
  }
}

getHeatMapCol <- function(mp, color.range = c("blue", "green", "yellow", "orange", "red")) {
  breaks <- round(mp, 0)
  if (all(is.na(breaks))) return(NULL)
  pal <- colorRampPalette(bias = 1, alpha = TRUE, colors = color.range, interpolate = "linear")(
    max(breaks, na.rm = TRUE) - min(breaks, na.rm = TRUE) + 1
  )
  cols <- data.frame(breaks = seq(min(breaks, na.rm = TRUE), max(breaks, na.rm = TRUE), 1),
                     col = pal, stringsAsFactors = FALSE)
  list(cols = cols$col[match(breaks, cols$breaks)], pal = pal)
}

color.bar <- function(lut, min, max = -min, nticks = 3, ticks = seq(min, max, len = nticks),
                      title = "", sideaxt = 4) {
  scale <- (length(lut) - 1) / (max - min)
  plot(c(0, 10), c(min, max), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = "")
  title(main = title, adj = 0, cex.main = .8, line = .5, xpd = TRUE)
  axis(sideaxt, round(ticks, 2), las = 1, line = -.9, tick = FALSE, cex.axis = .8, xpd = TRUE)
  for (i in 1:(length(lut) - 1)) {
    y <- (i - 1) / scale + min
    rect(0, y, 10, y + 1 / scale, col = lut[i], border = NA)
  }
}

addColorScale <- function(info, heatmap, colorTerm, xl = NULL, yl = NULL, sideaxt = 4, title = NULL) {
  xrange <- (par("usr")[2] - par("usr")[1]); yrange <- (par("usr")[4] - par("usr")[3])
  if (is.null(xl)) xl <- c(par("usr")[2] - xrange / 10, par("usr")[2] - xrange / 20)
  if (is.null(yl)) yl <- c(par("usr")[3] + yrange / 4, par("usr")[4] - yrange / 4)
  if (is.null(title)) title <- colorTerm
  Hmisc::subplot(
    fun = color.bar(lut = heatmap$pal,
                    min = min(info[, colorTerm], na.rm = TRUE),
                    title = title,
                    max = max(info[, colorTerm], na.rm = TRUE),
                    sideaxt = sideaxt),
    hadj = 1, vadj = 0, x = xl, y = yl
  )
}

MDS <- function(dist, info, group, groupbg, heat = FALSE, k = 2, main = NULL,
                ptSize = .6, pdfFile = NULL, cex.main = 1, cex.lab = 1, cex.axis = 1,
                mar = c(5, 4, 4, 2) + 0.1) {
  
  # distances: file path (TASSEL IBS) or square matrix/data.frame
  if (is.character(dist) && file.exists(dist)) {
    distances <- as.data.frame(data.table::fread(dist, sep = "\t", skip = 5, stringsAsFactors = FALSE))
    rownames(distances) <- distances$V1
    distances <- distances[, 2:ncol(distances)]
    colnames(distances) <- rownames(distances)
    if (nrow(distances) != ncol(distances)) R.utils::printf("Cannot correctly read %s! Expect TASSEL IBS.", dist)
  } else if (nrow(dist) == ncol(dist)) {
    distances <- dist
  } else {
    R.utils::printf("\nNon-valid entry for dist: must be TASSEL IBS file or square matrix/data.frame")
    return(NULL)
  }
  
  # info: file path or data.frame with "Taxon"
  if (is.character(info) && file.exists(info)) col <- read.table(info, header = TRUE, sep = "\t", as.is = TRUE) else col <- info
  if ("taxa" %in% tolower(names(col))) names(col)[tolower(names(col)) == "taxa"] <- "Taxon"
  if (!("taxon" %in% tolower(names(col)))) { R.utils::printf("\nInfo must have a 'Taxon' column!"); return(NULL) }
  taxonCol <- which("taxon" == tolower(names(col)))
  
  if (!(group %in% names(col))) { R.utils::printf("\nGroup variable not found in info."); return(NULL) }
  if (nrow(col) < 1) { R.utils::printf("\nNo rows in info file!"); return(NULL) }
  
  # if heat=FALSE and "<group>_color" exists, drop rows with empty color
  if (!heat) {
    colname <- paste0(group, "_color")
    if (colname %in% names(col)) col <- col[!is.na(col[[colname]]) & col[[colname]] != "", ]
  } else {
    col <- col[!(is.na(col[[group]]) | col[[group]] == ""), ]
  }
  
  # align distances & info
  distances <- distances[rownames(distances) %in% col[, taxonCol], rownames(distances) %in% col[, taxonCol]]
  if (nrow(distances) < 1) { R.utils::printf("\nNo rows left in distances after filtering."); return(NULL) }
  col <- col[col[, taxonCol] %in% rownames(distances), ]
  
  # iteratively drop samples driving NA columns
  rem.na <- apply(distances, 2, function(x) sum(is.na(x)))
  keep <- which(rem.na == 0 | rem.na < max(rem.na))
  while (length(keep) < nrow(distances)) {
    if ((nrow(distances) - length(keep)) > 0) {
      R.utils::printf("\nRemoving samples with NAs: %s",
                      paste(rownames(distances)[!(seq_len(nrow(distances)) %in% keep)], collapse = ","))
      distances <- distances[keep, keep]
    }
    rem.na <- apply(distances, 2, function(x) sum(is.na(x)))
    keep <- which(rem.na == 0 | rem.na < max(rem.na))
  }
  col <- col[col[, taxonCol] %in% rownames(distances), ]
  
  # metric MDS
  d <- as.matrix(distances)
  fit <- cmdscale(d, k = k, eig = TRUE)
  
  # colors
  if (heat) {
    heatmap <- getHeatMapCol(col[, group], color.range = c("blue", "yellow", "red"))
    color.acc <- heatmap$cols[match(rownames(fit$points[]), col[, taxonCol])]
  } else {
    legend <- levels(factor(col[, group][match(rownames(fit$points[]), col[, taxonCol])], levels = unique(col[, group])))
    if (!(paste0(group, "_color") %in% names(col))) {
      color.legend <- uniqueCols(length(legend))
      color.acc <- color.legend[match(col[, group], legend)]
    } else {
      color.legend <- col[, paste0(group, "_color")][match(legend, col[, group])]
      color.acc <- col[, paste0(group, "_color")][match(rownames(fit$points[]), col[, taxonCol])]
    }
  }
  
  # optional quick plot (kept as in source)
  if (!is.null(pdfFile)) svg(pdfFile, width = 6, height = 5)
  for (i in 1:(k - 1)) {
    par(mar = mar); j <- i + 1
    x <- fit$points[, i]; y <- fit$points[, j]
    if (is.null(main)) main <- paste("Metric MDS by", group)
    plot(x, y,
         xlab = paste0("Coordinate ", i, " (", round((fit$eig[i] / sum(fit$eig)) * 100, 1), "%)"),
         ylab = paste0("Coordinate ", j, " (", round((fit$eig[j] / sum(fit$eig)) * 100, 1), "%)"),
         main = main, col = color.acc, bg = color.acc, pch = 21, lwd = 2, cex = ptSize,
         xlim = c(min(fit$points[, 1]),
                  max(fit$points[, 1]) + (max(fit$points[, 1]) - min(fit$points[, 1])) * .3),
         cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)
    if (heat) {
      xrange <- (par("usr")[2] - par("usr")[1]); yrange <- (par("usr")[4] - par("usr")[3])
      addColorScale(info = col, heatmap = heatmap, colorTerm = group,
                    xl = c(par("usr")[2] - xrange / 5, par("usr")[2] - xrange / 6),
                    yl = c(par("usr")[3] + yrange / 4, par("usr")[4] - yrange / 4))
    } else {
      legend("bottomright", legend = legend, col = color.legend, pch = 19, cex = 0.6)
    }
  }
  if (!is.null(pdfFile)) dev.off()
  
  list(fit = fit, info = col, d = distances)
}

# ============================ Data Loading ====================================

load_data <- function() {
  # Distance matrix: TASSEL IBS format (5 header lines)
  distances <- as.data.frame(fread(distance_matrix, sep = "\t", skip = 5, stringsAsFactors = FALSE))
  rownames(distances) <- distances$V1
  distances <- distances[, 2:ncol(distances)]
  colnames(distances) <- rownames(distances)
  
  # Key file: supports comma decimals via dec=","
  dataMa <- as.data.frame(fread(MDS_KEY, sep = "\t", dec = ","))
  
  # Optional quality filter on 'proportion missing'
  dataMa$`proportion missing` <- as.numeric(dataMa$`proportion missing`)
  dataMa <- dataMa[dataMa$`proportion missing` < propmis, ]
  
  list(allDist = distances, dataMa = dataMa)
}

# compute MDS + return coordinates with standardized column names
create_mds_plot <- function(allDist, dataMa, cov_cutoff = NULL) {  # cov_cutoff kept for API compatibility
  mds_result <- MDS(dist = allDist, info = dataMa, group = color_group)
  mds_coords <- as.data.frame(mds_result$fit$points)
  colnames(mds_coords)[1:2] <- c("X1", "X2")
  mds_coords$Sample <- rownames(mds_coords)
  list(mds_result = mds_result, mds_coords = mds_coords)
}

# ============================ Shiny UI ========================================

ui <- fluidPage(
  titlePanel("Interactive plot for extracting sample names using polygons"),
  leafletOutput("mdsPlot")
)

# ============================ Shiny Server ====================================

server <- function(input, output, session) {
  
  # Load data & compute MDS once
  data <- load_data()
  mds_data <- create_mds_plot(data$allDist, data$dataMa)
  mds_coords <- mds_data$mds_coords
  
  # Build sf object for selection; CRS kept as 4326 to match your current approach
  mds_coords_sf <- st_as_sf(mds_coords, coords = c("X1", "X2"), crs = 4326)
  
  # Base map
  output$mdsPlot <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      addCircleMarkers(data = mds_coords_sf, radius = 1, color = "red", fillOpacity = 0.8) %>%
      addDrawToolbar(polyline = FALSE, polygon = TRUE)
  })
  
  # Group file counter (reactive)
  polygon_counter <- reactiveVal(1)
  
  # Handle drawn polygons
  observeEvent(input$mdsPlot_draw_new_feature, ignoreInit = TRUE, {
    drawn <- input$mdsPlot_draw_new_feature
    if (is.null(drawn)) return()
    
    # Build polygon in same CRS
    coords <- matrix(unlist(drawn$geometry$coordinates[[1]]), ncol = 2, byrow = TRUE)
    drawn_poly <- st_sf(geometry = st_sfc(st_polygon(list(coords))), crs = 4326)
    
    # Select points
    sel <- sf::st_intersection(mds_coords_sf, drawn_poly)
    
    # Highlight selection
    leafletProxy("mdsPlot") %>%
      clearMarkers() %>%
      addCircleMarkers(data = sel, radius = 5, color = "blue", fillOpacity = 0.8) %>%
      addCircleMarkers(data = mds_coords_sf, radius = 1, color = "red", fillOpacity = 0.8)
    
    # Export selected IDs
    ids <- unique(sel$Sample)
    message(sprintf("Selected sample names: %s", paste(ids, collapse = ", ")))
    
    if (!dir.exists(wd)) dir.create(wd, recursive = TRUE)
    out_file <- file.path(wd, sprintf("group%d.txt", polygon_counter()))
    writeLines(ids, out_file)
    polygon_counter(polygon_counter() + 1)
  })
}

# ============================ Run App =========================================

# sf warning about planar assumption is expected (MDS coordinates are not geographic).
shinyApp(ui, server)
