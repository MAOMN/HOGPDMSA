library(ComplexHeatmap)
library(circlize)
library(grid)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_2/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

# z-scores, groups and test results of the 13 mitochondrial genes (01_MitochondrialGenes_Analysis.R)
plotData <- readRDS( paste0( objectPath, "MitochondrialGenes_PlotData.rds" ) )
z_expr <- plotData$z_expr
group <- plotData$group
mtRes <- plotData$mtRes


#
### Heatmap (Figure 2A): z-scores of the samples and log2 fold change next to it
#

logfc <- matrix( mtRes$log2FoldChange, ncol = 1, dimnames = list( mtRes$symbol, "log2FC" ) )

# sizes in points
cellWidth = 11
cellHeight = 11.5
gapSize = unit( 2.85, "pt" )
barWidth = 8.5
barHeight = 61

zLimit <- max( 2, min( 2.5, max( abs( z_expr ) ) ) )
lfcLimit <- max( 0.5, ceiling( max( abs( logfc ) ) * 10 ) / 10 )

z_color <- colorRamp2( c( -zLimit, 0, zLimit ), c( "#6FA8DC", "white", "#E06666" ) )
lfc_color <- colorRamp2( c( -lfcLimit, 0, lfcLimit ), c( "#B39BE8", "white", "#F2CC5E" ) )

# ComplexHeatmap draws the colour bar 2 pt beyond the first and the last break (filled with the end colours),
# the bar is cut at the breaks: the ends of the bar are the limits of the colour scale and carry the end labels, as in the ggplot colour bars
barAtBreaks <- function( lg )
{
    body = getGrob( lg@grob, childNames( lg@grob )[ 2 ] )
    parts = childNames( body )
    colours = getGrob( body, parts[ 2 ] )
    ticks = convertY( getGrob( body, parts[ 3 ] )$y0, "pt", valueOnly = TRUE )
    n = length( colours$y )
    ends = c( 1, n )
    # the two end pieces lie outside the outer ticks
    stopifnot( length( parts ) == 5, abs( convertY( colours$y[ 1 ], "pt", valueOnly = TRUE ) - convertHeight( colours$height[ 1 ], "pt", valueOnly = TRUE ) / 2 - max( ticks ) ) < 1e-6,
               abs( convertY( colours$y[ n ], "pt", valueOnly = TRUE ) + convertHeight( colours$height[ n ], "pt", valueOnly = TRUE ) / 2 - min( ticks ) ) < 1e-6 )
    colours = editGrob( colours, x = colours$x[ -ends ], y = colours$y[ -ends ], height = colours$height[ -ends ], gp = gpar( col = colours$gp$col[ -ends ], fill = colours$gp$fill[ -ends ] ) )
    border = editGrob( getGrob( body, parts[ 5 ] ), y = unit( mean( range( ticks ) ), "pt" ), height = unit( diff( range( ticks ) ), "pt" ) )
    body = setGrob( setGrob( body, parts[ 2 ], colours ), parts[ 5 ], border )
    stopifnot( length( colours$x ) == n - 2, length( colours$y ) == n - 2, length( colours$gp$fill ) == n - 2 )
    lg@grob = setGrob( lg@grob, childNames( lg@grob )[ 2 ], body )
    lg
}

# color bars below each other
colorBar <- function( colFun, title, at ) barAtBreaks( Legend( col_fun = colFun, title = title, at = at, labels = cleanNumber( at ),
                                                  title_gp = gpar( fontfamily = "Arial", fontsize = axisTitleSize ),
                                                  labels_gp = gpar( fontfamily = "Arial", fontsize = legendValueSize ),
                                                  legend_height = unit( barHeight, "pt" ), grid_width = unit( barWidth, "pt" ),
                                                  border = "black", title_gap = unit( 6, "pt" ) ) )
legends <- packLegend( colorBar( z_color, "Z-score", seq( -zLimit, zLimit, length.out = 5 ) ),
                       colorBar( lfc_color, expression( log[2] * FC ), seq( -lfcLimit, lfcLimit, length.out = 5 ) ),
                       direction = "vertical", gap = unit( 6, "pt" ) )

# the groups are named above the columns
groupLabel <- factor( group, levels = levels( group ), labels = c( "Control", "aSyn" ) )

# the groups are named above the columns
groupLabel <- factor( group, levels = levels( group ), labels = c( "Control", "aSyn" ) )

hm_z <- Heatmap( z_expr, name = "Z-score", col = z_color,
                 width = ncol( z_expr ) * unit( cellWidth, "pt" ), height = nrow( z_expr ) * unit( cellHeight, "pt" ),
                 cluster_columns = FALSE, column_split = groupLabel, column_gap = gapSize,
                 column_title_gp = gpar( fontfamily = "Arial", fontsize = titleSize ), show_column_names = FALSE,
                 cluster_rows = FALSE, show_row_dend = FALSE, show_row_names = TRUE, row_names_side = "left",
                 row_names_gp = gpar( fontfamily = "Arial", fontface = "italic", fontsize = labelSize ),
                 rect_gp = gpar( col = NA ), border = TRUE, border_gp = gpar( col = "black", lwd = lineWidth ), show_heatmap_legend = FALSE )

# significance stars are written into the fold change column
hm_lfc <- Heatmap( logfc, name = "log2FC", col = lfc_color,
                   width = unit( cellWidth, "pt" ), height = nrow( logfc ) * unit( cellHeight, "pt" ),
                   cluster_columns = FALSE, show_column_names = FALSE, cluster_rows = FALSE, show_row_dend = FALSE, show_row_names = FALSE,
                   rect_gp = gpar( col = NA ), border = TRUE, border_gp = gpar( col = "black", lwd = lineWidth ),
                   cell_fun = function( j, i, x, y, width, height, fill ) {
                     if ( nzchar( mtRes$star[ i ] ) ) grid.text( mtRes$star[ i ], x, y, gp = gpar( fontfamily = "Arial", fontsize = annotationSize, col = "black" ) ) }, show_heatmap_legend = FALSE )

ht_opt$HEATMAP_LEGEND_PADDING <- unit( 8, "pt" )
ht_opt$TITLE_PADDING <- unit( 3, "pt" )

# title of the panel above the heatmaps, left aligned
title <- textGrob( "Mitochondrial protein-coding genes", x = 0, y = unit( 1, "npc" ) - unit( 2, "pt" ), hjust = 0, vjust = 1, gp = gpar( fontfamily = "Arial", fontsize = titleSize ) )
titleHeight <- unit( 1.2 * titleSize + 2, "pt" )

drawHeatmap <- function( newpage = TRUE ) draw( hm_z + hm_lfc, heatmap_legend_list = list( legends ), heatmap_legend_side = "right", ht_gap = gapSize,
                                                padding = unit( c( 4, 2, 2, 2 ), "pt" ), newpage = newpage )

# size of the heatmaps with legends, then of the whole figure
htSize <- local( {
    f = tempfile( fileext = ".pdf" )
    openPdf( f, 10, 10 )
    on.exit( { dev.off(); unlink( f ) } )
    ht = drawHeatmap()
    c( convertWidth( ComplexHeatmap:::width( ht ), "in", valueOnly = TRUE ), convertHeight( ComplexHeatmap:::height( ht ), "in", valueOnly = TRUE ),
       convertWidth( grobWidth( title ), "in", valueOnly = TRUE ) + 2 / 72, convertHeight( titleHeight, "in", valueOnly = TRUE ) )
} )
size <- c( max( htSize[ 1 ], htSize[ 3 ] ), htSize[ 2 ] + htSize[ 4 ] )

drawFigure2A <- function()
{
    pushViewport( viewport( x = 0, y = 1, width = unit( size[ 1 ], "in" ), height = unit( htSize[ 4 ], "in" ), just = c( "left", "top" ) ) )
    grid.draw( editGrob( title, x = unit( 2, "pt" ) ) )
    popViewport()
    pushViewport( viewport( x = 0, y = 0, width = unit( htSize[ 1 ], "in" ), height = unit( htSize[ 2 ], "in" ), just = c( "left", "bottom" ) ) )
    drawHeatmap( newpage = FALSE )
    popViewport()
}

saveFigure( drawFigure2A, paste0( figurePath, "Figure_2A_Heatmap_MitochondrialGenes" ), width = size[ 1 ], height = size[ 2 ], dpi = 600 )
