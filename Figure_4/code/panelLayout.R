library(ggplot2)
library(scales)
library(grid)
library(gtable)

# sizes and helper functions for the single panels of the main figures:
# text in Arial regular, all lines 0.5 pt, the plotting areas have a fixed size in points
# saveFigure.R has to be sourced before this file, a copy of both files is in each folder that needs them

# text sizes by kind of text
titleSize = 10         # panel titles, facet titles, group names above a heatmap
axisTitleSize = 10     # axis titles and legend titles
labelSize = 10         # names at the axes and in the legends: genes, transcription factors, cell states, gene sets
tickTextSize = 9       # numbers at the axes
legendValueSize = 8    # numbers in the legends
annotationSize = 9     # statistics and counts written into the panels
barTextSize = 7        # names written into bars
lineWidth = 0.5        # axes, box, ticks, curves
tickLength = 2
stripFill = "#B5B5B6"  # grey boxes with the names of the comparisons above the panels

# ggplot2: line widths and point sizes are not given in points
lw = lineWidth / ( .pt * 72 / 96 )
diameterToSize <- function( d ) d / ( 0.75 * .pt )

# theme of the panels: box or open axes, all text regular
themePanel <- function( box = TRUE )
{
    theme_classic( base_family = "Arial", base_size = axisTitleSize ) %+replace%
    theme( text = element_text( family = "Arial", colour = "black" ),
           plot.title = element_text( size = titleSize, hjust = 0.5, margin = margin( b = 3 ) ),
           axis.title = element_text( size = axisTitleSize ), axis.text = element_text( colour = "black", size = tickTextSize ),
           axis.line = if( box ) element_blank() else element_line( colour = "black", linewidth = lw ),
           axis.ticks = element_line( colour = "black", linewidth = lw ), axis.ticks.length = unit( tickLength, "pt" ),
           panel.border = if( box ) element_rect( colour = "black", fill = NA, linewidth = lw ) else element_blank(),
           strip.background = element_blank(), strip.text = element_text( size = titleSize ),
           legend.title = element_text( size = axisTitleSize ), legend.text = element_text( size = legendValueSize ),
           legend.background = element_blank(), legend.key = element_blank(), legend.margin = margin( 0, 0, 0, 0 ), legend.location = "plot",
           panel.background = element_blank(), plot.background = element_rect( fill = "white", colour = NA ), plot.margin = margin( 4, 4, 4, 4 ) )
}

# numbers without trailing zeros after the decimal point, zero written as 0
cleanNumber <- function( x )
{
    res = format( x, trim = TRUE, scientific = FALSE, drop0trailing = TRUE )
    res[ is.finite( x ) & abs( x ) < sqrt( .Machine$double.eps ) ] = "0"
    res[ is.na( x ) ] = NA
    res
}

# text is measured on the pdf device, the sizes are those of the pdf
withPdfDevice <- function( expr )
{
    f = tempfile( fileext = ".pdf" )
    openPdf( f, 10, 10 )
    on.exit( { dev.off(); unlink( f ) } )
    expr
}

# the plotting areas get a size in points
setPanelSize <- function( p, width, height )
{
    gt = withPdfDevice( ggplotGrob( p ) )
    panels = gt$layout[ grepl( "^panel", gt$layout$name ), ]
    gt$widths[ unique( panels$l ) ] <- unit( width, "pt" )
    gt$heights[ unique( panels$t ) ] <- unit( height, "pt" )
    # the size is fixed, a fixed aspect ratio (coord_equal) would shift the layout
    gt$respect <- FALSE
    # a legend at the side that reaches below the plot: the figure gets higher at the bottom
    # (the numbers at the ends of a colour bar reach half a line over the bar)
    for( i in which( gt$layout$name %in% c( "guide-box-right", "guide-box-left" ) ) )
    {
        missing = withPdfDevice( convertHeight( grobHeight( gt$grobs[[ i ]] ), "pt", valueOnly = TRUE ) -
                                   convertHeight( sum( gt$heights[ gt$layout$t[ i ]:nrow( gt ) ] ), "pt", valueOnly = TRUE ) ) + legendValueSize / 2 + 2
        if( missing > 0 )
        {
            gt = gtable_add_rows( gt, unit( missing, "pt" ), pos = -1 )
            gt$layout$b[ i ] = nrow( gt )
        }
    }
    gt
}

# size of a gtable in inch
tableSize <- function( gt ) withPdfDevice( c( convertWidth( sum( gt$widths ), "in", valueOnly = TRUE ), convertHeight( sum( gt$heights ), "in", valueOnly = TRUE ) ) )

# save as pdf, png and svg, the files have the size of the figure
saveGrob <- function( gt, fileName, dpi = 600 )
{
    size = tableSize( gt )
    saveFigure( function() grid.draw( gt ), fileName, width = size[ 1 ], height = size[ 2 ], dpi = dpi )
}

# pdf with every element as vector graphic (the files with points as image are lighter to edit)
saveFullVector <- function( gt, fileName )
{
    size = tableSize( gt )
    openPdf( paste0( fileName, "_FullVector.pdf" ), size[ 1 ], size[ 2 ] )
    grid.draw( gt )
    dev.off()
}

# save a plot with very many points: only the points are stored as image, axes, text, labels and legend stay vector graphics
# the same plot with the points as vector graphic is saved as well
savePointsAsImage <- function( gt, fileName, dpi = 600 )
{
    saveFullVector( gt, fileName )
    size = tableSize( gt )
    isPanel = grepl( "^panel", gt$layout$name )
    isPointLayer <- function( g ) grepl( "^geom_point", names( g$children ) )
    keepChildren <- function( g, keep ) { g$children <- g$children[ keep ]; g$childrenOrder <- names( g$children ); g }

    # only the point layers of the panels
    gPoints = gt
    for( i in which( !isPanel ) ) gPoints$grobs[[ i ]] = zeroGrob()
    for( i in which( isPanel ) ) gPoints$grobs[[ i ]] = keepChildren( gPoints$grobs[[ i ]], isPointLayer( gPoints$grobs[[ i ]] ) )
    # everything else
    gRest = gt
    for( i in which( gt$layout$name == "background" ) ) gRest$grobs[[ i ]] = zeroGrob()
    for( i in which( isPanel ) ) gRest$grobs[[ i ]] = keepChildren( gRest$grobs[[ i ]], !isPointLayer( gRest$grobs[[ i ]] ) )

    tmpFile = tempfile( "points_" )
    openPdf( paste0( tmpFile, ".pdf" ), size[ 1 ], size[ 2 ] )
    grid.newpage()
    grid.draw( gPoints )
    dev.off()
    system2( "pdftocairo", c( "-png", "-singlefile", "-r", dpi, shQuote( paste0( tmpFile, ".pdf" ) ), shQuote( tmpFile ) ) )
    pointImage = png::readPNG( paste0( tmpFile, ".png" ) )
    unlink( paste0( tmpFile, c( ".pdf", ".png" ) ) )

    saveFigure( function() {
      grid.newpage()
      grid.raster( pointImage, x = 0.5, y = 0.5, width = unit( 1, "npc" ), height = unit( 1, "npc" ), interpolate = FALSE )
      grid.draw( gRest ) }, fileName, width = size[ 1 ], height = size[ 2 ], dpi = dpi, pngFromPdf = TRUE )
}
