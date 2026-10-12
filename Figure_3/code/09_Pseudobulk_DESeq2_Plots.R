library(data.table)
library(ggplot2)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_3/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
tablePath = paste0( outputPath, "tables/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

# number of differentially expressed genes of 08_Pseudobulk_DESeq2_Analysis.R (the numbers of genes are integers)
summaryTab = as.data.frame( fread( paste0( tablePath, "Pseudobulk_DESeq2/Pseudobulk_DESeq2_Overview.txt" ) ) )

states = c( "OPC", paste0( "OL", 1:6 ) )


#
### Number of differentially expressed protein-coding genes per cell state (Figure 3D)
#

plotTab = rbind( data.frame( natural_label = summaryTab$natural_label, contrast = summaryTab$contrast, direction = "Up",
                             n_genes = summaryTab$n_up_lfc03_protein_coding, signed_genes = summaryTab$n_up_lfc03_protein_coding ),
                 data.frame( natural_label = summaryTab$natural_label, contrast = summaryTab$contrast, direction = "Down",
                             n_genes = summaryTab$n_down_lfc03_protein_coding, signed_genes = -summaryTab$n_down_lfc03_protein_coding ) )
plotTab$natural_label <- factor( plotTab$natural_label, levels = states )
plotTab$direction <- factor( plotTab$direction, levels = c( "Up", "Down" ) )
plotTab$contrast_label <- factor( ifelse( plotTab$contrast == "MSA_vs_HC", "MSA vs HC", "PD vs HC" ), levels = c( "MSA vs HC", "PD vs HC" ) )
plotTab <- plotTab[ order( plotTab$direction, plotTab$contrast_label, plotTab$natural_label ), ]

# y-axis: up-regulated genes upwards, down-regulated genes downwards, ticks in steps of 500 (or less if there are only few genes)
upMax = max( plotTab$signed_genes, 0 )
downMax = abs( min( plotTab$signed_genes, 0 ) )
yLimit = max( upMax, downMax )
yStep = if( yLimit <= 50 ) 10 else if( yLimit <= 200 ) 50 else if( yLimit <= 800 ) 200 else 500
upLimit = max( yStep, ceiling( upMax / yStep ) * yStep )
downLimit = max( yStep, ceiling( downMax / yStep ) * yStep )

# number of genes above (below) the bars
labelOffset = max( 5, max( upLimit, downLimit ) * 0.026 )
plotTab$label_y <- ifelse( plotTab$direction == "Up", ifelse( plotTab$n_genes > 0, plotTab$signed_genes + labelOffset, labelOffset ),
                           ifelse( plotTab$n_genes > 0, plotTab$signed_genes - labelOffset, -labelOffset ) )
plotTab$label_vjust <- ifelse( plotTab$direction == "Up", 0, 1 )

# size of one plotting area in points
panelWidth = 125
panelHeight = 149.5

p <- ggplot( plotTab, aes( natural_label, signed_genes, fill = direction ) ) +
  geom_hline( yintercept = 0, linewidth = lw, colour = "black" ) +
  geom_col( width = 0.8, colour = "black", linewidth = lw ) +
  geom_text( aes( y = label_y, label = n_genes, vjust = label_vjust ), family = "Arial", size = annotationSize / .pt, colour = "black", show.legend = FALSE ) +
  facet_grid( . ~ contrast_label ) +
  scale_fill_manual( values = c( Up = "#E06666", Down = "#6FA8DC" ), breaks = c( "Up", "Down" ), name = NULL ) +
  scale_y_continuous( breaks = seq( -downLimit, upLimit, by = yStep ), labels = function( x ) cleanNumber( abs( x ) ), limits = c( -downLimit * 1.18, upLimit * 1.18 ), expand = expansion( mult = c( 0, 0 ) ) ) +
  scale_x_discrete( limits = states, drop = FALSE, expand = expansion( mult = c( 0.1, 0.1 ) ) ) +
  labs( x = NULL, y = "DEG counts" ) + coord_cartesian( clip = "off" ) +
  themePanel( box = TRUE ) +
  theme( strip.background = element_rect( fill = stripFill, colour = "black", linewidth = lw ), strip.text = element_text( size = titleSize, margin = margin( 2, 0, 2, 0 ) ),
         axis.text.x = element_text( size = labelSize, angle = 45, hjust = 1, vjust = 1 ), panel.spacing.x = unit( 6, "pt" ),
         legend.position = "right", legend.text = element_text( size = labelSize ), legend.key.height = unit( 12, "pt" ), legend.key.width = unit( 12, "pt" ),
         legend.key.spacing.y = unit( 1, "pt" ) )

g <- setPanelSize( p, panelWidth, panelHeight )
saveGrob( g, paste0( figurePath, "Figure_3D_NumberOfDEGs_PerCellState" ) )
