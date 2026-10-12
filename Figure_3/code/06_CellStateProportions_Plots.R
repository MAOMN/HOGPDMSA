library(data.table)
library(ggplot2)
library(scales)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_3/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
tablePath = paste0( outputPath, "tables/CellStateProportions/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

# mean proportions of the cell states and the test results of 05_CellStateProportions_Analysis.R
plotData <- readRDS( paste0( objectPath, "Propeller_PlotData.rds" ) )
composition <- plotData$composition
statistics <- plotData$statistics

diagnoses = c( "HC", "MSA", "PD" )


#
### Stacked bar plot of the mean proportions (Figure 3C)
#

# positions of the bars: one bar per diagnosis, states stacked from OPC to OL6
composition$ymin <- 0
composition$ymax <- 0
for( dx in diagnoses )
{
    cand = which( composition$diagnosis == dx )
    composition$ymax[ cand ] = cumsum( composition$mean_prop[ cand ] )
    composition$ymin[ cand ] = cumsum( c( 0, composition$mean_prop[ cand ] )[ 1:length( cand ) ] )
}
composition$ymid <- ( composition$ymin + composition$ymax ) / 2
composition$label <- percent( composition$mean_prop, accuracy = 1 )
# no label for very small segments
composition$label_visible <- composition$mean_prop >= 0.035
composition$x_num <- as.numeric( composition$diagnosis )
composition$xmin <- composition$x_num - 0.32
composition$xmax <- composition$x_num + 0.32
# white text on the dark colors
composition$text_tone <- ifelse( as.character( composition$final_cluster_label ) %in% c( "OL1", "OL4", "OL6" ), "light", "dark" )

# test results of the MSA and PD bars (in the table, the plot shows no stars)
cand = match( paste( composition$diagnosis, composition$final_cluster_label ), paste( statistics$diagnosis, statistics$state ) )
composition$comparison <- statistics$comparison[ cand ]
composition$P.Value <- statistics$P.Value[ cand ]
composition$FDR <- statistics$FDR[ cand ]
composition$significance <- statistics$significance[ cand ]
composition$significance[ composition$diagnosis == "HC" ] <- ""
composition$stat_x <- composition$xmax + 0.13
composition$stat_y <- composition$ymid

fwrite( composition, file = paste0( tablePath, "OPC_OL_Composition.txt" ), sep = "\t" )

state_color = c( OPC = "#F2CC5E", OL1 = "#6FA8DC", OL2 = "#72C5C9", OL3 = "#95D5B2", OL4 = "#B39BE8", OL5 = "#D9826B", OL6 = "#3B86C4" )

# size of the plotting area in points
panelWidth = 116.5
panelHeight = 159.4

p <- ggplot( composition ) +
  geom_rect( aes( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = final_cluster_label ), colour = "black", linewidth = lw ) +
  geom_text( data = composition[ composition$label_visible & composition$text_tone == "dark", ], aes( x = x_num, y = ymid, label = label ),
             family = "Arial", size = annotationSize / .pt, colour = "black", show.legend = FALSE ) +
  geom_text( data = composition[ composition$label_visible & composition$text_tone == "light", ], aes( x = x_num, y = ymid, label = label ),
             family = "Arial", size = annotationSize / .pt, colour = "white", show.legend = FALSE ) +
  scale_fill_manual( values = state_color, name = "Cluster", drop = FALSE ) +
  scale_y_continuous( labels = function( x ) cleanNumber( 100 * x ), limits = c( 0, 1 ), expand = expansion( mult = c( 0, 0.035 ) ) ) +
  scale_x_continuous( breaks = 1:3, labels = diagnoses, limits = c( 0.55, 3.55 ), expand = expansion( mult = c( 0, 0 ) ) ) +
  labs( x = NULL, y = "Mean donor proportion (%)" ) +
  themePanel( box = TRUE ) +
  theme( axis.text.x = element_text( size = labelSize ), axis.ticks.x = element_blank(),
         panel.grid.major.y = element_line( colour = "#EFEDE8", linewidth = lw ),
         legend.position = "right", legend.text = element_text( size = labelSize ),
         legend.key.height = unit( 12, "pt" ), legend.key.width = unit( 12, "pt" ), legend.key.spacing.y = unit( 1, "pt" ) )

g <- setPanelSize( p, panelWidth, panelHeight )
saveGrob( g, paste0( figurePath, "Figure_3C_CellStateProportions" ) )
