library(ggplot2)
library(scales)
library(grid)
library(gtable)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_2/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

# running enrichment scores, ranked Wald statistics and test results of 03_GSEA_MitoCarta_Analysis.R
plotData <- readRDS( paste0( objectPath, "GSEA_PlotData.rds" ) )

# sizes of the plotting areas in points: running score, gene positions, ranked statistics
panelWidth = 156
panelHeights = c( 80, 12, 34 )


# enrichment curve with gene positions and ranked statistics below
plotEnrichmentCurve <- function( curve, rankTab, resRow, title, outFile )
{
    hitTab = curve[ curve$hit, ]
    leadTab = curve[ curve$leading_edge, ]

    lineColor = if( resRow$NES >= 0 ) "#E06666" else "#6FA8DC"
    fdrText = if( resRow$padj >= 0.01 ) sprintf( "%.2f", resRow$padj ) else scientific_format( digits = 2 )( resRow$padj )
    statText = paste0( "NES = ", sprintf( "%.2f", resRow$NES ), "\n", "FDR = ", fdrText, "\n", "LE = ", nrow( leadTab ), "/", nrow( hitTab ) )
    # statistics in the upper right corner
    statPos = data.frame( rank = nrow( rankTab ) * 0.98, running_es = Inf, label = statText )

    xLimits = c( 1, nrow( rankTab ) )
    xBreaks = seq( 2000, xLimits[ 2 ], by = 4000 )
    stacked = theme( plot.margin = margin( 0, 4, 0, 4 ) )

    # running score
    p_es <- ggplot( curve, aes( rank, running_es ) ) + geom_hline( yintercept = 0, linewidth = lw, colour = "black" ) +
      geom_line( colour = lineColor, linewidth = lw ) +
      geom_text( data = statPos, aes( rank, running_es, label = label ), inherit.aes = FALSE, family = "Arial", size = annotationSize / .pt, lineheight = 0.92, hjust = 1, vjust = 1.15, colour = "black" ) +
      scale_x_continuous( limits = xLimits, breaks = xBreaks, expand = expansion( mult = c( 0, 0 ) ) ) +
      scale_y_continuous( labels = cleanNumber, expand = expansion( mult = c( 0.08, 0.11 ) ) ) +
      labs( x = NULL, y = "Enrichment score" ) + themePanel( box = TRUE ) + stacked +
      theme( axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin( 0, 4, 0, 4 ) )

    # positions of the genes of the gene set
    p_hits <- ggplot( hitTab, aes( rank, 0 ) ) +
      geom_segment( aes( xend = rank, y = 0, yend = 1, colour = leading_edge ), linewidth = lw ) +
      scale_colour_manual( values = c( `TRUE` = lineColor, `FALSE` = "black" ), guide = "none" ) +
      scale_x_continuous( limits = xLimits, breaks = xBreaks, expand = expansion( mult = c( 0, 0 ) ) ) +
      scale_y_continuous( limits = c( 0, 1 ), expand = expansion( mult = c( 0, 0 ) ) ) +
      labs( x = NULL, y = NULL ) + themePanel( box = TRUE ) + stacked +
      theme( axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank() )

    # ranked Wald statistics
    p_rank <- ggplot( rankTab, aes( rank, stat_plot ) ) + geom_hline( yintercept = 0, linewidth = lw, colour = "black" ) +
      geom_area( data = rankTab[ rankTab$stat_plot >= 0, ], aes( y = stat_plot ), fill = "#E06666", alpha = 0.78, linewidth = 0 ) +
      geom_area( data = rankTab[ rankTab$stat_plot < 0, ], aes( y = stat_plot ), fill = "#6FA8DC", alpha = 0.78, linewidth = 0 ) +
      scale_x_continuous( limits = xLimits, breaks = xBreaks, labels = cleanNumber, expand = expansion( mult = c( 0, 0 ) ) ) +
      scale_y_continuous( limits = range( rankTab$stat_plot ), breaks = pretty_breaks( n = 3 ), labels = cleanNumber ) +
      labs( x = "Ranked protein-coding genes", y = "Wald stat" ) + themePanel( box = TRUE ) + theme( plot.margin = margin( 0, 4, 4, 4 ) )

    # the three plots below each other without space, the title above
    parts = mapply( setPanelSize, list( p_es, p_hits, p_rank ), panelWidth, panelHeights, SIMPLIFY = FALSE )
    g = do.call( rbind, c( parts, size = "max" ) )
    titleGrob = textGrob( title, gp = gpar( fontfamily = "Arial", fontsize = titleSize, lineheight = 0.95 ) )
    g = gtable_add_rows( g, unit( 1, "grobheight", titleGrob ) + unit( 7, "pt" ), pos = 0 )
    panels = g$layout[ grepl( "^panel", g$layout$name ), ]
    g = gtable_add_grob( g, titleGrob, t = 1, l = min( panels$l ), r = max( panels$r ), clip = "off" )

    saveGrob( g, outFile )
}


# Figure 2B: all protein-coding genes, Figure 2C: without the genes of the mitochondrial genome
plotEnrichmentCurve( plotData$all$curve, plotData$all$rankPlot, plotData$all$resRow, "MitoCarta: OXPHOS subunits",
                     paste0( figurePath, "Figure_2B_GSEA_OXPHOS_Subunits" ) )
plotEnrichmentCurve( plotData$nuclear$curve, plotData$nuclear$rankPlot, plotData$nuclear$resRow, "Nuclear-coding genes MitoCarta:\nOXPHOS subunits",
                     paste0( figurePath, "Figure_2C_GSEA_OXPHOS_Subunits_NuclearGenes" ) )
