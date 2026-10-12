library(ggplot2)
library(ggrastr)
library(scales)
library(viridisLite)
library(gtable)

# change to your project directory
home_directory <- "path/to/project_directory/"

# output folder
outputPath = paste0( home_directory, "Figure_5/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )


# two genes of the respiratory chain that are in the first module of PD L1 and of PD L3: COX6B1 and ATP5F1E
disease = "PD"
genes = c( "COX6B1", "ATP5F1E" )

# fitted curves and cells of the two genes (made by 05_SelectedGenes_OXPHOS_PD_Analysis.R)
plotData <- readRDS( paste0( objectPath, "SelectedGenes_OXPHOS_", disease, "_PlotData.rds" ) )
curveTab = plotData$curveTab
pointTab = plotData$pointTab


#
### Plot (Figure 5D and E): one figure per lineage, both genes next to each other
#

# size of the plotting areas and of the points, in points
panelWidth = 108.4
panelHeight = 92.1
pointDiameter = 1.3
pointSize = pointDiameter / ( 0.75 * .pt )

condition_color = setNames( viridis( 2 ), c( "HC", "Disease" ) )

# y-axis starts a bit below zero so that the cells without expression are fully visible
yLimits <- function( limits )
{
    upper = max( limits, na.rm = TRUE )
    if( !is.finite( upper ) || upper <= 0 ) upper = 1
    c( -0.08 * upper, upper )
}
yBreaks <- function( limits )
{
    b = breaks_pretty( n = 5 )( c( 0, max( limits, na.rm = TRUE ) ) )
    b[ is.finite( b ) & b >= 0 ]
}
yLabels <- function( x )
{
    res = label_number( accuracy = 0.1, trim = TRUE )( x )
    res = sub( "[.]0+$", "", res )
    res = sub( "([.][0-9]*[1-9])0+$", "\\1", res )
    res[ abs( x ) < .Machine$double.eps^0.5 ] = "0"
    res
}

for( lin in c( "L1", "L3" ) )
{
    cv = curveTab[ curveTab$lineage_id == lin, ]
    pt = pointTab[ pointTab$lineage_id == lin, ]
    cv$gene = factor( cv$gene, levels = genes )
    pt$gene = factor( pt$gene, levels = genes )
    cv$condition_plot = factor( cv$condition_plot, levels = c( "HC", "Disease" ) )
    pt$condition_plot = factor( pt$condition_plot, levels = c( "HC", "Disease" ) )

    # the cells as image (and as vector graphic in the file *_FullVector.pdf)
    for( rasterized in c( TRUE, FALSE ) )
    {
        pointLayer = if( rasterized ) geom_point_rast( data = pt, aes( x = pseudotime, y = log_expression, colour = condition_plot ), size = pointSize, alpha = 0.75, stroke = 0, raster.dpi = 600, dev = "cairo" ) else
                                      geom_point( data = pt, aes( x = pseudotime, y = log_expression, colour = condition_plot ), size = pointSize, alpha = 0.75, stroke = 0 )
        p <- ggplot() + pointLayer +
          geom_line( data = cv, aes( x = pseudotime, y = fitted_log_expression, colour = condition_plot, group = condition_plot ), linewidth = lw, lineend = "round" ) +
          facet_wrap( vars( gene ), ncol = 2, scales = "free_y" ) +
          scale_y_continuous( limits = yLimits, breaks = yBreaks, labels = yLabels, expand = expansion( mult = c( 0, 0.04 ) ) ) +
          scale_colour_manual( values = condition_color, breaks = c( "HC", "Disease" ), labels = c( "HC", disease ), name = NULL ) +
          labs( title = paste0( disease, " vs HC ", lin ), x = "OL pseudotime", y = "Ln(expression + 1)" ) +
          guides( colour = guide_legend( override.aes = list( alpha = 1, linewidth = lw, size = pointSize ) ) ) +
          theme_classic( base_family = "Arial", base_size = axisTitleSize ) +
          theme( text = element_text( colour = "black" ),
                 plot.title = element_text( size = titleSize, hjust = 0.5, margin = margin( b = 3 ) ),
                 axis.title = element_text( size = axisTitleSize ), axis.title.y = element_text( margin = margin( r = 5 ) ), axis.text = element_text( colour = "black", size = tickTextSize ),
                 axis.line = element_line( colour = "black", linewidth = lw ), axis.ticks = element_line( colour = "black", linewidth = lw ), axis.ticks.length = unit( tickLength, "pt" ),
                 strip.background = element_blank(), strip.text = element_text( size = labelSize, face = "italic", margin = margin( t = 4, b = 2 ) ), panel.spacing = unit( 5, "pt" ),
                 legend.position = "bottom", legend.text = element_text( size = labelSize ), legend.key.width = unit( 15, "pt" ), legend.key.height = unit( 10, "pt" ),
                 legend.margin = margin( 0, 0, 0, 0 ), legend.box.spacing = unit( 3, "pt" ),
                 plot.margin = margin( 6, 4, 3, 4 ) )

        gt = setPanelSize( p, panelWidth, panelHeight )

        # rule below the title, from the left end of the first panel to the right end of the second
        panels = gt$layout[ grepl( "^panel", gt$layout$name ), ]
        titleRow = gt$layout$t[ gt$layout$name == "title" ]
        gt <- gtable_add_grob( gt, segmentsGrob( x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar( lwd = lineWidth * 96 / 72, col = "black", lineend = "butt" ) ),
                               t = titleRow, l = min( panels$l ), r = max( panels$r ), clip = "off", name = "title_rule" )

        fileName = paste0( figurePath, "Figure_5", if( lin == "L1" ) "D" else "E", "_SelectedGenes_OXPHOS_", disease, "_", lin )
        if( rasterized ) saveGrob( gt, fileName ) else saveFullVector( gt, fileName )
    }
}

# add the two panels to the list of the panel files of Figure 5 (made by 04_Disease_Dynamics_Modules_Plots.R, run it first)
panelFiles = data.frame( panel = c( "D", "E" ), disease = disease, lineage = c( "L1", "L3" ), module_id = c( "PD_L1__D01", "PD_L3__D01" ),
                         figure_pdf = c( "Figure_5D_SelectedGenes_OXPHOS_PD_L1.pdf", "Figure_5E_SelectedGenes_OXPHOS_PD_L3.pdf" ) )
write.table( panelFiles, file = paste0( figurePath, "Figure_5_PanelFiles.txt" ), append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
