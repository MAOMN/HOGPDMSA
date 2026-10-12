library(data.table)
library(ggplot2)
library(scales)
library(grid)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_3/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
tablePath = paste0( outputPath, "tables/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )


# AUCell scores of 10_AUCell_OXPHOS_Analysis.R
scores = as.data.frame( fread( paste0( tablePath, "AUCell/AUCell_OXPHOS_Signature_Scores.txt.gz" ) ) )
scores$diagnosis <- factor( scores$diagnosis, levels = c( "HC", "MSA", "PD" ) )

# cells with high score are drawn last (on top), cells with the same score in random order
set.seed( 42 )
scores <- scores[ order( frank( scores$AUC, ties.method = "random" ) ), ]

# color scale from 0 to the 99% quantile of the scores (rounded up), larger values get the darkest color
q99 = as.numeric( quantile( scores$AUC, 0.99 ) )
stepRaw = q99 / 3
magnitude = 10^floor( log10( stepRaw ) )
niceSteps = c( 1, 1.5, 2, 2.5, 5, 10 )
colorStep = niceSteps[ which( niceSteps >= stepRaw / magnitude )[ 1 ] ] * magnitude
colorMax = ceiling( q99 / colorStep - 1e-12 ) * colorStep
colorBreaks = seq( 0, colorMax, by = colorStep )
# number of decimals that are needed for the step
nDigits = which( sapply( 0:8, function( k ) abs( colorStep * 10^k - round( colorStep * 10^k ) ) < 1e-08 ) )[ 1 ] - 1
colorLabels = formatC( colorBreaks, format = "f", digits = nDigits )
# no zero at the end of a number with decimals (0.10 is written 0.1)
colorLabels = ifelse( grepl( "\\.", colorLabels ), sub( "\\.?0+$", "", colorLabels ), colorLabels )
colorLabels[ 1 ] = "0"
colorLabels[ length( colorLabels ) ] = paste0( "≥", colorLabels[ length( colorLabels ) ] )

# axis ticks in steps of 5 within the range of the cells
umapBreaks <- function( x )
{
    b = seq( 5 * ceiling( min( x ) / 5 ), 5 * floor( max( x ) / 5 ), by = 5 )
    b[ b >= min( x ) & b <= max( x ) ]
}
xBreaks = umapBreaks( scores$UMAP_1 )
yBreaks = umapBreaks( scores$UMAP_2 )
xLimits = range( scores$UMAP_1 )
yLimits = range( scores$UMAP_2 )

# sizes of the plot, in points
aucPanelWidth = 129.9
# the plotting area has the same scale on both axes
aucPanelHeight = aucPanelWidth * diff( yLimits ) / diff( xLimits )
aucPointDiameter = 1
aucBarWidth = 7.7
aucBarHeight = 52.8

plotAUC <- function( data, split = FALSE )
{
    p <- ggplot( data, aes( UMAP_1, UMAP_2, color = AUC ) ) +
      geom_point( size = diameterToSize( aucPointDiameter ), alpha = if( split ) 0.92 else 0.93, stroke = 0 ) +
      scale_color_gradientn( colours = c( "#3B86C4", "#D7E5EF", "#F6F4EF", "#E4A0A6", "#8F2638" ), limits = c( 0, colorMax ), breaks = colorBreaks,
                             labels = colorLabels, oob = squish, name = "AUC" ) +
      scale_x_continuous( breaks = xBreaks, limits = xLimits, expand = expansion( mult = c( 0.025, 0.025 ) ) ) +
      scale_y_continuous( breaks = yBreaks, limits = yLimits, expand = expansion( mult = c( 0.025, 0.025 ) ) ) +
      coord_equal( clip = "off" ) + labs( x = "UMAP 1", y = "UMAP 2" ) +
      themePanel( box = TRUE ) +
      theme( legend.position = "right", legend.direction = "vertical", strip.text = element_text( size = titleSize ) ) +
      guides( color = guide_colorbar( title.position = "top", title.hjust = 0, barheight = unit( aucBarHeight, "pt" ), barwidth = unit( aucBarWidth, "pt" ),
                                      frame.colour = "black", frame.linewidth = lw, ticks.colour = "black", ticks.linewidth = lw ) )

    if( split ) p <- p + facet_wrap( ~diagnosis, nrow = 1, drop = FALSE ) + theme( panel.spacing.x = unit( 11, "pt" ) )
    setPanelSize( p, aucPanelWidth, aucPanelHeight )
}

### Figure 3E: HC, MSA and PD next to each other
savePointsAsImage( plotAUC( scores, split = TRUE ), paste0( figurePath, "Figure_3E_AUCell_OXPHOS_Signature_HC_MSA_PD" ) )

### same as single plots
for( dx in c( "HC", "MSA", "PD" ) ) savePointsAsImage( plotAUC( scores[ scores$diagnosis == dx, ] ), paste0( figurePath, "AUCell_OXPHOS_Signature_", dx ) )
