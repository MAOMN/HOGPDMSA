library(ggplot2)
library(grid)
library(gtable)

# functions for the dot plots of transcription factor activities
# - one square per transcription factor and comparison
# - a dot is drawn if the activity is significant: color = VIPER NES, size = -log10 FDR


#
# Numbers in the legends without trailing zeros, zero is written as 0
#
legendNumber <- function( x )
{
    sapply( x, function( v ) if( is.na( v ) ) NA_character_ else if( v == 0 ) "0" else format( v, trim = TRUE, scientific = FALSE, digits = 12 ) )
}


#
# Round the largest value up to a nice upper limit for the legend (1, 1.5, 2, ..., 5, 10 times a power of ten)
#
niceLegendLimit <- function( x )
{
    x = abs( x[ is.finite( x ) ] )
    if( length( x ) == 0 || max( x ) == 0 )
    {
        return( 1 )
    }
    top = max( x )
    decade = 10^floor( log10( top ) )
    cand = decade * c( 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 10 )
    return( min( cand[ cand >= top ] ) )
}


#
# Breaks of the size legend
#
legendSizeBreaks <- function( top )
{
    decade = 10^floor( log10( top ) )
    step = if( top / decade <= 5 ) 0.5 * decade else 5 * decade
    res = seq( step, top, by = step )
    if( length( res ) > 5 )
    {
        res = res[ unique( round( seq( 1, length( res ), length.out = 5 ) ) ) ]
    }
    return( res )
}


theme_TF <- function()
{
    theme_minimal( base_family = "Arial", base_size = 11 ) +
      theme( panel.grid = element_blank(), panel.border = element_rect( fill = NA, colour = "black", linewidth = 0.5 ),
             axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_text( colour = "black" ),
             axis.text.y = element_text( face = "bold.italic", size = 11, margin = margin( r = 4 ) ),
             axis.text.x = element_text( face = "bold", size = 10.5, angle = 45, hjust = 1, vjust = 1 ),
             strip.background = element_rect( fill = "#F5F5F5", colour = "black", linewidth = 0.4 ), strip.text = element_text( face = "bold", size = 11 ),
             legend.title = element_text( face = "bold", size = 10.5 ), legend.text = element_text( size = 10 ),
             legend.box = "vertical", legend.position = "bottom", legend.direction = "horizontal", legend.spacing.y = unit( 1, "mm" ),
             legend.box.spacing = unit( 2, "mm" ), panel.spacing.x = unit( 7, "pt" ), plot.margin = margin( 2, 5, 2, 5 ) )
}


#
# Size and color scale of the dots
# scoreMax, fdrMax - upper limits of the legends
# narrow - size legend in two rows (for small plots)
# fiveTicks - 5 instead of 3 numbers at the color bar
#
scales_TF <- function( scoreMax, fdrMax, narrow = FALSE, fiveTicks = FALSE )
{
    sizeBreaks = legendSizeBreaks( fdrMax )
    colorBreaks = if( fiveTicks ) c( -scoreMax, -scoreMax / 2, 0, scoreMax / 2, scoreMax ) else c( -scoreMax, 0, scoreMax )
    list( scale_size_continuous( name = expression( bold( paste( "−log"[10], " FDR" ) ) ), range = c( 1.6, 3.6 ), limits = c( 0, fdrMax ), breaks = sizeBreaks, labels = legendNumber ),
          scale_fill_gradient2( name = "VIPER NES", low = "#6FA8DC", mid = "#F7F4EF", high = "#E06666", midpoint = 0, limits = c( -scoreMax, scoreMax ), breaks = colorBreaks, labels = legendNumber ),
          guides( size = guide_legend( order = 1, title.position = "top", nrow = if( narrow && length( sizeBreaks ) > 3 ) 2 else 1, byrow = TRUE, override.aes = list( fill = "#C7C7C7", shape = 21 ) ),
                  fill = guide_colourbar( order = 2, title.position = "top", barwidth = unit( 28, "mm" ), barheight = unit( 3, "mm" ), frame.colour = "black", frame.linewidth = 0.5,
                                          ticks = TRUE, ticks.colour = "black", ticks.linewidth = 0.5 ) ) )
}


#
# All squares get the same size (0.215 inch), the size of the figure follows from the number of rows and columns
# returns the plot as grob together with its width and height in inch
#
fixedSquares <- function( p, nRows, nCols )
{
    g = ggplotGrob( p )
    panels = g$layout[ grepl( "^panel", g$layout$name ), ]
    g$widths[ unique( panels$l ) ] = unit( nCols * 0.215, "in" )
    g$heights[ unique( panels$t ) ] = unit( nRows * 0.215, "in" )
    w = max( convertWidth( sum( g$widths ), "in", valueOnly = TRUE ), 2.10 )
    h = convertHeight( sum( g$heights ), "in", valueOnly = TRUE )
    return( list( g = g, w = w, h = h ) )
}
