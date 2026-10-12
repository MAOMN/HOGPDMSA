library(data.table)
library(ggplot2)
library(grid)
library(gtable)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_3/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
tablePath = paste0( outputPath, "tables/TF_Activity/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "plotTFActivity.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

# activities and ranking of the transcription factors of 14_TF_Activity_OPC_OL_Analysis.R
tfTable = as.data.frame( fread( paste0( tablePath, "TF_Activity_OPC_OL_AllTFs.txt" ) ) )
ranking = as.data.frame( fread( paste0( tablePath, "TF_Activity_OPC_OL_Ranking_PD_FDRSum.txt" ) ) )

states = c( "OPC", paste0( "OL", 1:6 ) )


#
### Plot (Figure 3G): 10 transcription factors that are most consistently significant over the cell states in PD
#

topTF = ranking$source[ 1:10 ]

plotTab <- tfTable[ tfTable$source %in% topTF, ]
plotTab$significant <- !is.na( plotTab$FDR ) & plotTab$FDR < 0.05
plotTab$neglog10_FDR <- ifelse( plotTab$significant, -log10( pmax( plotTab$FDR, .Machine$double.xmin ) ), NA )
plotTab$source_plot <- factor( plotTab$source, levels = rev( topTF ) )
plotTab$stratum_plot <- factor( plotTab$stratum, levels = states )
plotTab$contrast_plot <- factor( plotTab$disease, levels = c( "MSA", "PD" ), labels = c( "MSA vs HC", "PD vs HC" ) )

sigPoints = plotTab[ plotTab$significant, ]
scoreMax = niceLegendLimit( abs( sigPoints$score ) )
fdrMax = niceLegendLimit( sigPoints$neglog10_FDR )

# sizes of the dot plot, in points
tfCell = 12.8                      # one square
tfDotRange = tfCell * c( 0.125, 0.53 ) # diameter of the dots at the lower and the upper limit of the size scale
tfDotStroke = 0.5 * 96 / 72 / ( .stroke / 2 )
tfBarWidth = 64
tfBarHeight = 6.8

p <- ggplot( plotTab, aes( stratum_plot, source_plot ) ) +
  geom_tile( width = 1, height = 1, fill = "white", colour = "black", linewidth = lw ) +
  geom_point( data = sigPoints, aes( size = neglog10_FDR, fill = score ), shape = 21, colour = "black", stroke = tfDotStroke ) +
  facet_wrap( vars( contrast_plot ), nrow = 1, drop = FALSE ) +
  scale_x_discrete( expand = expansion( add = 0 ) ) + scale_y_discrete( expand = expansion( add = 0 ) ) +
  scale_size_continuous( name = expression( paste( "−log"[10], " FDR" ) ), range = diameterToSize( tfDotRange ), limits = c( 0, fdrMax ),
                         breaks = legendSizeBreaks( fdrMax ), labels = legendNumber ) +
  scale_fill_gradient2( name = "TF activity", low = "#6FA8DC", mid = "#F7F4EF", high = "#E06666", midpoint = 0, limits = c( -scoreMax, scoreMax ),
                        breaks = c( -scoreMax, -scoreMax / 2, 0, scoreMax / 2, scoreMax ), labels = legendNumber ) +
  guides( size = guide_legend( order = 2, title.position = "top", nrow = 1, byrow = TRUE, override.aes = list( fill = "#C7C7C7", shape = 21 ) ),
          fill = guide_colourbar( order = 1, title.position = "top", barwidth = unit( tfBarWidth, "pt" ), barheight = unit( tfBarHeight, "pt" ), frame.colour = "black", frame.linewidth = lw,
                                  ticks = TRUE, ticks.colour = "black", ticks.linewidth = lw ) ) +
  themePanel( box = TRUE ) +
  theme( axis.title = element_blank(), axis.ticks = element_blank(),
         axis.text.y = element_text( face = "italic", size = labelSize, colour = "black", margin = margin( r = 2 ) ),
         axis.text.x = element_text( size = labelSize, colour = "black", angle = 45, hjust = 1, vjust = 1 ),
         strip.background = element_rect( fill = stripFill, colour = "black", linewidth = lw ), strip.text = element_text( size = titleSize, margin = margin( 2, 0, 2, 0 ) ),
         legend.title = element_text( size = axisTitleSize ), legend.text = element_text( size = legendValueSize ),
         legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", legend.box.spacing = unit( 4, "pt" ),
         legend.spacing.x = unit( 14, "pt" ), legend.key.size = unit( 10, "pt" ), panel.spacing.x = unit( 5.8, "pt" ) )

g <- setPanelSize( p, length( states ) * tfCell, length( topTF ) * tfCell )

saveGrob( g, paste0( figurePath, "Figure_3G_TF_Activity_OPC_OL" ) )
