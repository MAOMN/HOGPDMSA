library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(GO.db)

# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/" )

# output folder
outputPath = paste0( home_directory, "Figure_4/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

lineageTab = data.frame( unit_id = c( "HC_L1_TDE", "HC_L2_TDE", "HC_L3_TDE" ), trajectory_id = c( "L1", "L2", "L3" ), model_lineage = 1:3, endpoint = c( "OL5", "OL6", "OL3" ) )

# fitted expression curves of the genes that change along pseudotime in HC (one file per lineage, made by 08_HC_Dynamics_AssociationTest_Analysis.R)
# list with: genes, raw_curve (genes x 50 bins), row_z (z-scores), tested_genes, bins
curves = list()
for( i in 1:3 )
{
    curves[[ lineageTab$unit_id[ i ] ]] = readRDS( paste0( objectPath, "HC_Dynamics_Curves_", lineageTab$trajectory_id[ i ], ".rds" ) )
}

# modules, GO terms and distribution of the cell clusters along pseudotime (made by 09_HC_Dynamics_Modules_Analysis.R)
# list with: topTab (genes shown in the heatmaps), goTop3 (GO terms written next to the heatmaps), waves (cell clusters along pseudotime)
plotData <- readRDS( paste0( objectPath, "HC_Dynamics_PlotData.rds" ) )
topTab = plotData$topTab
goTop3 = plotData$goTop3
waves = plotData$waves


#
### Heatmaps (Figure 4B-D): one figure per lineage
#

# marker genes of the oligodendrocyte lineage that are labelled if they are in the heatmap
# columns: gene, stage, stage_order, display_order
markerGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/Oligodendrocyte_LineageMarkerGenes.txt" ) ) )

module_colors = c( "#C05A8B", "#7A8A1C", "#315A93", "#2EA6C7" )
heat_colors = c( "#315A93", "#3B86C4", "#6FA8DC", "#95D5B2", "#F2CC5E" )
cluster_color = c( OPC = "#F2CC5E", OL1 = "#6FA8DC", OL2 = "#E06666", OL3 = "#95D5B2", OL4 = "#B39BE8", OL5 = "#E58F73", OL6 = "#3B86C4" )

# sizes in points; two of the figures fit next to each other on one page (width of 612 pt)
heatWidth = 75        # heatmap
heatHeight = 185
waveHeight = 12       # one cell cluster above the heatmap
blockWidth = 5        # coloured bar of the modules
barWidth = 50         # colour bar of the z-scores
barHeight = 6.8

# names of the GO terms: official name of GO.db, first letter upper case
sentenceCase <- function( x ) paste0( toupper( substr( x, 1, 1 ) ), substring( x, 2 ) )
goName <- function( id, fallback ) ifelse( is.na( id ), fallback, sentenceCase( unname( Term( id ) ) ) )

textWidth <- function( x, size, face = "plain" ) withPdfDevice( max( sapply( x, function( v ) convertWidth( grobWidth( textGrob( v, gp = gpar( fontfamily = "Arial", fontsize = size, fontface = face ) ) ), "pt", valueOnly = TRUE ) ) ) )

#
# Positions of the marker gene labels left of the heatmap: labels need a minimal distance to each other
#
spreadLabels <- function( rowIndex, n )
{
    lower = max( 8, n * 0.04 )
    upper = min( n - 8, n * 0.96 )
    k = length( rowIndex )
    if( k == 1 )
    {
        return( pmin( pmax( rowIndex, lower ), upper ) )
    }
    separation = min( n * 0.055, ( upper - lower ) / ( k - 1 ) )
    index = 1:k
    projected = isoreg( index, rowIndex - ( index - 1 ) * separation )$yf
    projected = pmin( pmax( projected, lower ), upper - ( k - 1 ) * separation )
    return( projected + ( index - 1 ) * separation )
}

for( i in 1:3 )
{
    u = lineageTab$unit_id[ i ]
    top = topTab[ topTab$unit_id == u, ]
    mat = curves[[ u ]]$row_z[ top$gene, ]
    n = nrow( mat )

    #rows of each module
    ranges = c()
    for( b in unique( top$block ) )
    {
        rows = top$display_row_order[ top$block == b ]
        goRows = goTop3[ goTop3$unit_id == u & goTop3$block == b, ]
        go = goName( goRows$go_id, goRows$term_label )
        ranges = rbind( ranges, data.frame( block = b, block_order = top$block_order[ top$block == b ][ 1 ], ymin = min( rows ) - 0.5, ymax = max( rows ) + 0.5, ymid = mean( range( rows ) ),
                                            displayed_genes = length( rows ),
                                            go_label = paste( sapply( go, function( x ) paste( strwrap( x, width = 32 ), collapse = "\n" ) ), collapse = "\n" ) ) )
    }
    ranges$colour = module_colors[ ranges$block_order ]

    #cell clusters along pseudotime
    wave = waves[[ u ]]$smooth
    wave$final_cluster_label = factor( wave$final_cluster_label, levels = waves[[ u ]]$path )
    nClusters = length( waves[[ u ]]$path )
    p_wave <- ggplot( wave, aes( pseudotime_scaled, count_scaled, fill = final_cluster_label ) ) +
      geom_area( colour = "black", linewidth = lw, alpha = 0.96 ) + facet_grid( final_cluster_label ~ . ) +
      scale_fill_manual( values = cluster_color, guide = "none", drop = FALSE ) +
      scale_x_continuous( limits = c( 0, 1 ), breaks = NULL, expand = expansion( mult = c( 0, 0 ) ) ) +
      scale_y_continuous( limits = c( 0, 1.08 ), breaks = NULL, expand = expansion( mult = c( 0, 0 ) ) ) +
      labs( x = NULL, y = NULL, title = paste( "HC Lineage", sub( "^L", "", lineageTab$trajectory_id[ i ] ) ) ) + themePanel( box = TRUE ) +
      theme( axis.ticks = element_blank(), strip.text = element_blank(), panel.spacing.y = unit( 0, "pt" ), plot.margin = margin( 0, 0, 3, 0 ) )

    # names of the cell clusters left of their rows (in the column of the marker genes)
    p_waveName <- ggplot( data.frame( final_cluster_label = factor( waves[[ u ]]$path, levels = waves[[ u ]]$path ) ), aes( 1, 0.5, label = final_cluster_label ) ) +
      geom_text( family = "Arial", size = labelSize / .pt, hjust = 1, nudge_x = -0.03 ) + facet_grid( final_cluster_label ~ . ) +
      scale_x_continuous( limits = c( 0, 1 ), expand = expansion( mult = c( 0, 0 ) ) ) + scale_y_continuous( limits = c( 0, 1 ), expand = expansion( mult = c( 0, 0 ) ) ) +
      coord_cartesian( clip = "off" ) + theme_void( base_family = "Arial" ) +
      theme( strip.text = element_blank(), panel.spacing.y = unit( 0, "pt" ), plot.margin = margin( 0, 0, 3, 0 ) )

    #heatmap of the z-scores
    heat = data.frame( fit_position = rep( 1:ncol( mat ), times = n ), row_index = rep( top$display_row_order, each = ncol( mat ) ), zscore = as.vector( t( mat ) ) )
    p_heat <- ggplot( heat, aes( fit_position, row_index, fill = zscore ) ) + geom_raster() +
      geom_rect( data = ranges, aes( xmin = 0.5, xmax = ncol( mat ) + 0.5, ymin = ymin, ymax = ymax ), inherit.aes = FALSE, fill = NA, colour = "black", linewidth = lw ) +
      scale_fill_gradientn( colours = heat_colors, values = rescale( c( -2, -1.15, -0.25, 0.75, 2 ) ), limits = c( -2, 2 ), breaks = c( -2, 0, 2 ), labels = c( "-2", "0", "2" ), oob = squish, name = "Z-score" ) +
      scale_x_continuous( limits = c( 0.5, ncol( mat ) + 0.5 ), breaks = c( 0.5, ( ncol( mat ) + 1 ) / 2, ncol( mat ) + 0.5 ), labels = c( "0", "0.5", "1" ), expand = expansion( mult = c( 0, 0 ) ) ) +
      scale_y_reverse( limits = c( n + 0.6, 0.4 ), breaks = NULL, expand = expansion( mult = c( 0, 0 ) ) ) +
      labs( x = "Pseudotime", y = NULL ) + themePanel( box = TRUE ) +
      theme( axis.ticks.y = element_blank(), legend.position = "bottom", legend.title.position = "top", plot.margin = margin( 0, 0, 0, 0 ) ) +
      guides( fill = guide_colourbar( title.hjust = 0.5, barwidth = unit( barWidth, "pt" ), barheight = unit( barHeight, "pt" ),
                                      frame.colour = "black", frame.linewidth = lw, ticks.colour = "black", ticks.linewidth = lw ) )

    #marker genes
    markers = markerGenes[ markerGenes$gene %in% top$gene, ]
    markers$row_index = top$display_row_order[ match( markers$gene, top$gene ) ]
    markers = markers[ order( markers$stage_order, markers$display_order, markers$row_index ), ]
    markers = markers[ order( markers$row_index ), ]
    markers$label_y = spreadLabels( markers$row_index, n )
    markerWidth = textWidth( markers$gene, labelSize, "italic" ) + 16
    p_marker <- ggplot( markers ) +
      geom_segment( aes( x = 0.73, xend = 0.84, y = label_y, yend = label_y ), colour = "black", linewidth = lw ) +
      geom_segment( aes( x = 0.84, xend = 0.84, y = label_y, yend = row_index ), colour = "black", linewidth = lw ) +
      geom_segment( aes( x = 0.84, xend = 1, y = row_index, yend = row_index ), colour = "black", linewidth = lw ) +
      geom_text( aes( x = 0.7, y = label_y, label = gene ), family = "Arial", fontface = "italic", size = labelSize / .pt, hjust = 1, colour = "black" ) +
      scale_x_continuous( limits = c( 0, 1 ), expand = expansion( mult = c( 0, 0 ) ) ) + scale_y_reverse( limits = c( n + 0.6, 0.4 ), expand = expansion( mult = c( 0, 0 ) ) ) +
      coord_cartesian( clip = "off" ) + theme_void( base_family = "Arial" ) + theme( plot.margin = margin( 0, 0, 0, 0 ) )

    #colored bar of the modules
    p_block <- ggplot( ranges, aes( xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = block ) ) + geom_rect( colour = "black", linewidth = lw ) +
      scale_fill_manual( values = setNames( module_colors, ranges$block ), guide = "none" ) +
      scale_x_continuous( limits = c( 0, 1 ), expand = expansion( mult = c( 0, 0 ) ) ) + scale_y_reverse( limits = c( n + 0.6, 0.4 ), expand = expansion( mult = c( 0, 0 ) ) ) +
      theme_void( base_family = "Arial" ) + theme( plot.margin = margin( 0, 0, 0, 0 ) )

    #GO terms of the modules, in the colour of the module
    goWidth = textWidth( unlist( strsplit( ranges$go_label, "\n" ) ), labelSize ) + 6
    p_go <- ggplot( ranges ) +
      geom_text( aes( x = 0, y = ymid, label = go_label, colour = block ), family = "Arial", size = labelSize / .pt, hjust = 0, vjust = 0.5, lineheight = 0.9 ) +
      scale_colour_manual( values = setNames( module_colors, ranges$block ), guide = "none" ) +
      scale_x_continuous( limits = c( 0, 1 ), expand = expansion( add = c( 0.02, 0 ) ) ) + scale_y_reverse( limits = c( n + 0.6, 0.4 ), expand = expansion( mult = c( 0, 0 ) ) ) +
      coord_cartesian( clip = "off" ) + theme_void( base_family = "Arial" ) + theme( plot.margin = margin( 0, 0, 0, 0 ) )

    #"GO BP" above the GO terms
    p_goTitle <- ggplot() + annotate( "text", x = 0, y = 0, label = "GO BP", family = "Arial", size = titleSize / .pt, hjust = 0, vjust = 0 ) +
      scale_x_continuous( limits = c( 0, 1 ), expand = expansion( mult = c( 0, 0 ) ) ) + scale_y_continuous( limits = c( 0, 1 ), expand = expansion( mult = c( 0, 0 ) ) ) +
      coord_cartesian( clip = "off" ) + theme_void( base_family = "Arial" ) + theme( plot.margin = margin( 0, 0, 3, 3 ) )

    # one grid of 2 x 4 plots: the plotting areas of the columns and rows are aligned
    empty <- plot_spacer() + theme( plot.margin = margin( 0, 0, 0, 0 ) )
    p <- wrap_plots( p_waveName, p_wave, empty, p_goTitle, p_marker, p_heat, p_block, p_go, ncol = 4,
                     widths = unit( c( markerWidth, heatWidth, blockWidth, goWidth ), "pt" ), heights = unit( c( nClusters * waveHeight, heatHeight ), "pt" ) ) +
      plot_annotation( theme = theme( plot.margin = margin( 4, 4, 4, 4 ) ) )

    gt <- withPdfDevice( patchworkGrob( p ) )
    saveGrob( gt, paste0( figurePath, "Figure_4_HC_Dynamics_Heatmap_", lineageTab$trajectory_id[ i ] ) )
}
