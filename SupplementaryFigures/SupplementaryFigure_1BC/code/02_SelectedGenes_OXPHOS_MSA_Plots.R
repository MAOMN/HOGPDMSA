library(ggplot2)
library(ggrastr)
library(scales)
library(viridisLite)
library(grid)
library(gtable)

# change to your project directory
home_directory <- "path/to/project_directory/"

# output folder
outputPath = paste0( home_directory, "SupplementaryFigures/SupplementaryFigure_1BC/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )


# the two genes of the respiratory chain that are shown for PD in Figure 5D and E (COX6B1 and ATP5F1E), here for MSA
disease = "MSA"
genes = c( "COX6B1", "ATP5F1E" )

# fitted curves and cells of the two genes (made by 01_SelectedGenes_OXPHOS_MSA_Analysis.R)
plotData <- readRDS( paste0( objectPath, "SelectedGenes_OXPHOS_", disease, "_PlotData.rds" ) )
curveTab = plotData$curveTab
pointTab = plotData$pointTab

# sizes of the figure, all in points (the files have the final size of the figure)
textSize = 10          # axis titles, panel titles, legend
tickTextSize = 9       # numbers at the axes
letterSize = 13        # panel letters, bold
lineWidth = 0.5        # axes, ticks, fitted curves, rule below the header
tickLength = 2
panelWidth = 108.4
panelHeight = 92.1
pointDiameter = 1.3
# ggplot: linewidth and point size are not given in points
lw = lineWidth / ( .pt * 72 / 96 )
pointSize = pointDiameter / ( 0.75 * .pt )


#
### Plot (Supplementary Figure 1B and 1C): lineage 1 (B) and lineage 3 (C) next to each other, both genes of a lineage next to each other
#

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

#
# Both genes of one lineage next to each other, header with a rule below it and the panel letter
#
lineageBlock <- function( lin, letter, rasterized = TRUE )
{
    cv = curveTab[ curveTab$lineage_id == lin, ]
    pt = pointTab[ pointTab$lineage_id == lin, ]
    cv$gene = factor( cv$gene, levels = genes )
    pt$gene = factor( pt$gene, levels = genes )
    cv$condition_plot = factor( cv$condition_plot, levels = c( "HC", "Disease" ) )
    pt$condition_plot = factor( pt$condition_plot, levels = c( "HC", "Disease" ) )

    # the cells as image or as vector graphic
    pointLayer = if( rasterized ) geom_point_rast( data = pt, aes( x = pseudotime, y = log_expression, colour = condition_plot ), size = pointSize, alpha = 0.75, stroke = 0, raster.dpi = 600, dev = "cairo" ) else
                                  geom_point( data = pt, aes( x = pseudotime, y = log_expression, colour = condition_plot ), size = pointSize, alpha = 0.75, stroke = 0 )
    p <- ggplot() + pointLayer +
      geom_line( data = cv, aes( x = pseudotime, y = fitted_log_expression, colour = condition_plot, group = condition_plot ), linewidth = lw, lineend = "round" ) +
      facet_wrap( vars( gene ), ncol = 2, scales = "free_y" ) +
      scale_y_continuous( limits = yLimits, breaks = yBreaks, labels = yLabels, expand = expansion( mult = c( 0, 0.04 ) ) ) +
      scale_colour_manual( values = condition_color, breaks = c( "HC", "Disease" ), labels = c( "HC", disease ), name = NULL ) +
      labs( title = paste0( disease, " vs HC ", lin ), x = "OL pseudotime", y = "Ln(expression + 1)" ) +
      guides( colour = guide_legend( override.aes = list( alpha = 1, linewidth = lw, size = pointSize ) ) ) +
      theme_classic( base_family = "Arial", base_size = textSize ) +
      theme( text = element_text( colour = "black" ),
             plot.title = element_text( size = textSize, hjust = 0.5, margin = margin( b = 3 ) ),
             axis.title = element_text( size = textSize ), axis.title.y = element_text( margin = margin( r = 5 ) ), axis.text = element_text( colour = "black", size = tickTextSize ),
             axis.line = element_line( colour = "black", linewidth = lw ), axis.ticks = element_line( colour = "black", linewidth = lw ), axis.ticks.length = unit( tickLength, "pt" ),
             strip.background = element_blank(), strip.text = element_text( size = textSize, face = "italic", margin = margin( t = 4, b = 2 ) ), panel.spacing = unit( 5, "pt" ),
             legend.position = "bottom", legend.text = element_text( size = textSize ), legend.key.width = unit( 15, "pt" ), legend.key.height = unit( 10, "pt" ),
             legend.margin = margin( 0, 0, 0, 0 ), legend.box.spacing = unit( 3, "pt" ),
             plot.margin = margin( 6, 4, 3, 4 ) )

    gt = ggplotGrob( p )

    # all panels have the same size
    panels = gt$layout[ grepl( "^panel", gt$layout$name ), ]
    gt$widths[ unique( panels$l ) ] <- unit( panelWidth, "pt" )
    gt$heights[ unique( panels$t ) ] <- unit( panelHeight, "pt" )

    # rule below the header, from the left end of the first panel to the right end of the second
    titleRow = gt$layout$t[ gt$layout$name == "title" ]
    gt <- gtable_add_grob( gt, segmentsGrob( x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar( lwd = lineWidth * 96 / 72, col = "black", lineend = "butt" ) ),
                           t = titleRow, l = min( panels$l ), r = max( panels$r ), clip = "off", name = "header_rule" )

    # panel letter at the upper left corner (not for the single panels)
    if( !is.null( letter ) ) gt <- gtable_add_grob( gt, textGrob( letter, x = 0, y = 1, hjust = 0, vjust = 1, gp = gpar( fontfamily = "Arial", fontface = "bold", fontsize = letterSize, col = "black" ) ),
                           t = 2, b = titleRow, l = 2, r = min( panels$l ), clip = "off", name = "panel_letter" )
    gt
}

#
# Size of a gtable in inch
#
tableSize <- function( gt ) c( convertWidth( sum( gt$widths ), "in", valueOnly = TRUE ), convertHeight( sum( gt$heights ), "in", valueOnly = TRUE ) )

# Supplementary Figure 1B and 1C
blocks = list( lineageBlock( "L1", "B" ), lineageBlock( "L3", "C" ) )
gap = unit( 12, "pt" )
fig = gtable( widths = unit.c( sum( blocks[[ 1 ]]$widths ), gap, sum( blocks[[ 2 ]]$widths ) ), heights = max( sum( blocks[[ 1 ]]$heights ), sum( blocks[[ 2 ]]$heights ) ) )
fig <- gtable_add_grob( fig, blocks, t = 1, l = c( 1, 3 ) )
size = tableSize( fig )
saveFigure( function() grid.draw( fig ), paste0( figurePath, "SupplementaryFigure_1BC" ), width = size[ 1 ], height = size[ 2 ], dpi = 600 )

# single panels
single = lineageBlock( "L1", NULL )
size = tableSize( single )
saveFigure( function() grid.draw( single ), paste0( figurePath, "Panels/SupplementaryFigure_1B" ), width = size[ 1 ], height = size[ 2 ], dpi = 600 )
single = lineageBlock( "L3", NULL )
size = tableSize( single )
saveFigure( function() grid.draw( single ), paste0( figurePath, "Panels/SupplementaryFigure_1C" ), width = size[ 1 ], height = size[ 2 ], dpi = 600 )

# the same figures with the cells as vector graphic (pdf only)
saveFullVector <- function( g, fileName )
{
    size = tableSize( g )
    openPdf( paste0( fileName, "_FullVector.pdf" ), size[ 1 ], size[ 2 ] )
    grid.draw( g )
    dev.off()
}
blocks = list( lineageBlock( "L1", "B", rasterized = FALSE ), lineageBlock( "L3", "C", rasterized = FALSE ) )
fig = gtable( widths = unit.c( sum( blocks[[ 1 ]]$widths ), gap, sum( blocks[[ 2 ]]$widths ) ), heights = max( sum( blocks[[ 1 ]]$heights ), sum( blocks[[ 2 ]]$heights ) ) )
fig <- gtable_add_grob( fig, blocks, t = 1, l = c( 1, 3 ) )
saveFullVector( fig, paste0( figurePath, "SupplementaryFigure_1BC" ) )
saveFullVector( lineageBlock( "L1", NULL, rasterized = FALSE ), paste0( figurePath, "Panels/SupplementaryFigure_1B" ) )
saveFullVector( lineageBlock( "L3", NULL, rasterized = FALSE ), paste0( figurePath, "Panels/SupplementaryFigure_1C" ) )
