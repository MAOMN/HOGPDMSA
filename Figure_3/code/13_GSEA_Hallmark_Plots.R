library(stringr)
library(ggplot2)
library(scales)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_3/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

# GSEA results and the most significant gene sets of 12_GSEA_Hallmark_Analysis.R
plotData <- readRDS( paste0( objectPath, "GSEA_Hallmark_PlotData.rds" ) )
gseaTab <- plotData$gseaTab
topTerms <- plotData$topTerms

states = c( "OPC", paste0( "OL", 1:6 ) )


#
### Plot (Figure 3F): the 4 most significant gene sets of each cell state and comparison
#

plotTerms = unique( topTerms$term_id )

# order of the gene sets in the plot: clustering by signed significance (NES * -log10 FDR, 0 if not significant)
columns = sort( unique( paste( gseaTab$contrast, gseaTab$natural_label, sep = "_" ) ), method = "radix" )
scoreMat = matrix( 0, nrow = length( plotTerms ), ncol = length( columns ), dimnames = list( plotTerms, columns ) )
for( i in which( gseaTab$significant & gseaTab$term_id %in% plotTerms ) )
{
    scoreMat[ gseaTab$term_id[ i ], paste( gseaTab$contrast[ i ], gseaTab$natural_label[ i ], sep = "_" ) ] = gseaTab$NES[ i ] * gseaTab$neg_log10_fdr[ i ]
}
rowOrder = rownames( scoreMat )[ hclust( dist( scoreMat ), method = "ward.D2" )$order ]

termTab = data.frame( term_id = plotTerms, display_term = gseaTab$display_term[ match( plotTerms, gseaTab$term_id ) ] )
termTab$row_rank = match( termTab$term_id, rowOrder )
termTab = termTab[ order( termTab$row_rank ), ]
termTab$y = nrow( termTab ) - 1:nrow( termTab ) + 1
termTab$display_term_wrapped = str_wrap( termTab$display_term, width = 36 )

# x positions: MSA vs HC left, PD vs HC right with a gap in between
gapWidth = 0.7
columnTab = rbind( data.frame( contrast = "MSA_vs_HC", natural_label = states, x = 1:7 ),
                   data.frame( contrast = "PD_vs_HC", natural_label = states, x = 1:7 + ( 7 + gapWidth ) ) )

# one square per gene set, cell state and comparison
plotTab = c()
for( i in 1:nrow( termTab ) )
{
    d = columnTab[ order( columnTab$contrast, columnTab$natural_label ), ]
    d$term_id = termTab$term_id[ i ]
    d$y = termTab$y[ i ]
    cand = match( paste( d$term_id, d$contrast, d$natural_label ), paste( gseaTab$term_id, gseaTab$contrast, gseaTab$natural_label ) )
    d$NES = gseaTab$NES[ cand ]
    d$padj = gseaTab$padj[ cand ]
    d$significant = gseaTab$significant[ cand ]
    d$neg_log10_fdr = gseaTab$neg_log10_fdr[ cand ]
    plotTab = rbind( plotTab, d )
}
plotPoints = plotTab[ which( plotTab$significant ), ]

nTerms = nrow( termTab )
left_xmin = min( columnTab$x[ columnTab$contrast == "MSA_vs_HC" ] ) - 0.5
left_xmax = max( columnTab$x[ columnTab$contrast == "MSA_vs_HC" ] ) + 0.5
right_xmin = min( columnTab$x[ columnTab$contrast == "PD_vs_HC" ] ) - 0.5
right_xmax = max( columnTab$x[ columnTab$contrast == "PD_vs_HC" ] ) + 0.5
strip_ymin = nTerms + 0.6
strip_ymax = nTerms + 1.55

# color scale: symmetric, upper limit rounded up to the next half step
nesMax = max( abs( plotPoints$NES ) )
halfStep = 0.5 * 10^floor( log10( nesMax ) )
nesLimit = ceiling( nesMax / halfStep ) * halfStep
if( nesLimit - nesMax <= sqrt( .Machine$double.eps ) * max( 1, nesMax ) ) nesLimit = nesLimit + halfStep
nesBreaks = seq( -nesLimit, nesLimit, length.out = 5 )
nesDigits = max( 0, 1 - floor( log10( nesMax ) ) )
nesLabels = sapply( nesBreaks, function( v ) if( abs( v ) < sqrt( .Machine$double.eps ) ) "0" else formatC( v, format = "f", digits = nesDigits ) )
# no zero at the end of a number with decimals (-3.0 is written -3)
nesLabels = ifelse( grepl( "\\.", nesLabels ), sub( "\\.?0+$", "", nesLabels ), nesLabels )

sizeBreaks = c( 2, 5, 10, 20 )
sizeBreaks = sizeBreaks[ sizeBreaks >= floor( min( plotPoints$neg_log10_fdr ) ) & sizeBreaks <= ceiling( max( plotPoints$neg_log10_fdr ) ) ]

# one square in points, the plot has the same scale on both axes
cellSize = 11.1
maxDotDiameter = 0.9 * cellSize
barWidth = 64
barHeight = 6.8

p <- ggplot( plotTab, aes( x = x, y = y ) ) +
  geom_tile( width = 1, height = 1, fill = NA, colour = "black", linewidth = lw ) +
  geom_point( data = plotPoints, aes( size = neg_log10_fdr, fill = NES ), shape = 21, colour = "black", stroke = 0.5 * 96 / 72 / ( .stroke / 2 ), alpha = 0.96 ) +
  annotate( "rect", xmin = left_xmin, xmax = left_xmax, ymin = 0.5, ymax = nTerms + 0.5, fill = NA, colour = "black", linewidth = lw ) +
  annotate( "rect", xmin = right_xmin, xmax = right_xmax, ymin = 0.5, ymax = nTerms + 0.5, fill = NA, colour = "black", linewidth = lw ) +
  annotate( "rect", xmin = left_xmin, xmax = left_xmax, ymin = strip_ymin, ymax = strip_ymax, fill = stripFill, colour = "black", linewidth = lw ) +
  annotate( "rect", xmin = right_xmin, xmax = right_xmax, ymin = strip_ymin, ymax = strip_ymax, fill = stripFill, colour = "black", linewidth = lw ) +
  annotate( "text", x = ( left_xmin + left_xmax ) / 2, y = ( strip_ymin + strip_ymax ) / 2, label = "MSA vs HC", family = "Arial", size = titleSize / .pt, colour = "black" ) +
  annotate( "text", x = ( right_xmin + right_xmax ) / 2, y = ( strip_ymin + strip_ymax ) / 2, label = "PD vs HC", family = "Arial", size = titleSize / .pt, colour = "black" ) +
  scale_x_continuous( breaks = columnTab$x, labels = columnTab$natural_label, limits = c( left_xmin, right_xmax ), expand = expansion( mult = c( 0, 0 ) ) ) +
  scale_y_continuous( breaks = termTab$y, labels = termTab$display_term_wrapped, limits = c( 0.5, nTerms + 1.6 ), expand = expansion( mult = c( 0, 0 ) ) ) +
  scale_fill_gradient2( name = "NES", low = "#6FA8DC", mid = "#F7F4EF", high = "#E06666", midpoint = 0, limits = c( -nesLimit, nesLimit ), breaks = nesBreaks, labels = nesLabels,
                        guide = guide_colourbar( title.position = "top", barwidth = unit( barWidth, "pt" ), barheight = unit( barHeight, "pt" ), direction = "horizontal",
                                                 frame.colour = "black", frame.linewidth = lw, ticks.colour = "black", ticks.linewidth = lw, order = 1 ) ) +
  scale_size_area( name = expression( paste( "−log"[10], " FDR" ) ), max_size = diameterToSize( maxDotDiameter ), breaks = sizeBreaks, labels = number_format( accuracy = 1 ),
                   guide = guide_legend( title.position = "top", override.aes = list( fill = "#C7C7C7", colour = "black", alpha = 1 ), order = 2 ) ) +
  labs( x = NULL, y = NULL, title = NULL ) + coord_fixed( ratio = 1, clip = "off" ) +
  themePanel( box = FALSE ) +
  theme( axis.text.x = element_text( size = labelSize, angle = 45, hjust = 1, vjust = 1 ),
         axis.text.y = element_text( size = labelSize, lineheight = 0.9 ),
         axis.ticks = element_blank(), axis.line = element_blank(),
         legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = unit( 14, "pt" ), legend.key.size = unit( 10, "pt" ) )

g <- setPanelSize( p, ( right_xmax - left_xmin ) * cellSize, ( nTerms + 1.1 ) * cellSize )
saveGrob( g, paste0( figurePath, "Figure_3F_GSEA_Hallmark" ) )
