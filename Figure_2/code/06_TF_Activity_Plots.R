library(ggplot2)
library(grid)
library(gtable)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_2/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "plotTFActivity.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

# significant transcription factors of 05_TF_Activity_Analysis.R
sigTF <- readRDS( paste0( objectPath, "TF_Activity_SignificantTFs.rds" ) )


#
### Plot of the significant transcription factors (Figure 2F)
#

plotTab <- sigTF[ order( -sigTF$score, sigTF$source ), ]
plotTab$source_plot <- factor( plotTab$source, levels = plotTab$source )
plotTab$neglog10_FDR <- -log10( pmax( plotTab$FDR, .Machine$double.xmin ) )

scoreMax = niceLegendLimit( abs( plotTab$score ) )
fdrMax = niceLegendLimit( plotTab$neglog10_FDR )

# sizes of the dot plot, in points
tfCell = 15.5                          # one square
tfDotRange = tfCell * c( 0.125, 0.53 ) # diameter of the dots at the lower and the upper limit of the size scale
tfDotStroke = 0.5 * 96 / 72 / ( .stroke / 2 )
tfBarWidth = 64
tfBarHeight = 6.8

# the names at 45 degrees reach over the right end of the row, the numbers at the ends of the color bar over its ends
textWidth <- function( x, size, face = "plain" ) withPdfDevice( max( sapply( x, function( v ) convertWidth( grobWidth( textGrob( v, gp = gpar( fontfamily = "Arial", fontsize = size, fontface = face ) ) ), "pt", valueOnly = TRUE ) ) ) )
rightMargin = max( 4, textWidth( levels( plotTab$source_plot )[ nrow( plotTab ) ], labelSize, "italic" ) * cos( pi / 4 ) - tfCell / 2 + 4 )
sideMargin = 4 + textWidth( legendNumber( -scoreMax ), legendValueSize ) / 2

p <- ggplot( plotTab, aes( source_plot, 1 ) ) +
  geom_tile( width = 1, height = 1, fill = "white", colour = "black", linewidth = lw ) +
  geom_point( aes( size = neglog10_FDR, fill = score ), shape = 21, colour = "black", stroke = tfDotStroke ) +
  scale_x_discrete( position = "top", expand = expansion( add = 0 ) ) +
  scale_y_continuous( breaks = NULL, limits = c( 0.5, 1.5 ), expand = expansion( add = 0 ) ) +
  scale_size_continuous( name = expression( paste( "−log"[10], " FDR" ) ), range = diameterToSize( tfDotRange ), limits = c( 0, fdrMax ),
                         breaks = legendSizeBreaks( fdrMax ), labels = legendNumber ) +
  scale_fill_gradient2( name = "TF activity", low = "#6FA8DC", mid = "#F7F4EF", high = "#E06666", midpoint = 0, limits = c( -scoreMax, scoreMax ),
                        breaks = c( -scoreMax, -scoreMax / 2, 0, scoreMax / 2, scoreMax ), labels = legendNumber ) +
  guides( size = guide_legend( order = 2, title.position = "top", nrow = 1, byrow = TRUE, override.aes = list( fill = "#C7C7C7", shape = 21 ) ),
          fill = guide_colourbar( order = 1, title.position = "top", barwidth = unit( tfBarWidth, "pt" ), barheight = unit( tfBarHeight, "pt" ), frame.colour = "black", frame.linewidth = lw,
                                  ticks = TRUE, ticks.colour = "black", ticks.linewidth = lw ) ) +
  themePanel( box = TRUE ) +
  theme( axis.title = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
         axis.text.x.top = element_text( face = "italic", size = labelSize, colour = "black", angle = 45, hjust = 0, vjust = 0 ),
         legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", legend.box.spacing = unit( 8, "pt" ),
         legend.spacing.x = unit( 14, "pt" ), legend.key.size = unit( 10, "pt" ),
         plot.margin = margin( 4, max( rightMargin, sideMargin ), 4, sideMargin ) )

g <- setPanelSize( p, nrow( plotTab ) * tfCell, tfCell )

saveGrob( g, paste0( figurePath, "Figure_2F_TF_Activity_HOG" ) )
