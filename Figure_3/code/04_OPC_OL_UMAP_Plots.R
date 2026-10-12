library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_3/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
tablePath = paste0( outputPath, "tables/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )


# table of the clustered OPCs and oligodendrocytes (03_OPC_OL_Clustering_Analysis.R)
# columns: cell_id, UMAP_1, UMAP_2, nCount_RNA, nFeature_RNA, sample_id, diagnosis, labels_major, labels_minor, natural_cluster, final_cluster_label,
#          final_retained, small_cluster, donor_support_pass, n_donors_with_min_cells, final_exclusion_reason
cells = as.data.frame( fread( paste0( tablePath, "OPC_OL_Clustering/OPC_OL_CellTable.txt.gz" ) ) )

# cells of the kept clusters
keptCells <- cells[ cells$final_retained, ]

# colors
cluster_colors = c( "#F2CC5E", "#6FA8DC", "#E06666", "#95D5B2", "#B39BE8", "#D9826B", "#3B86C4", "#5DA271", "#C97AA8", "#7A5CCF",
                    "#B8A038", "#6C8EBF", "#A65E3A", "#2AA198", "#8E6C3A", "#C65D7B", "#4B5563", "#00A6A6", "#A64AC9", "#7AC36A" )
diagnosis_color = c( HC = "#D0CDD6", MSA = "#E06666", PD = "#6FA8DC" )

# width of the plotting area in points, the height follows from the same scale on both axes
umapWidth = 143.5
legendDotDiameter = 6


# UMAP colored by one column of the cell table
# labelOnData: write the names into the plot instead of a legend, otherwise the legend is in one row at the bottom of the plotting area
plotUMAP <- function( data, colorBy, palette, labelOnData = FALSE )
{
    palette = palette[ names( palette ) %in% unique( as.character( data[ , colorBy ] ) ) ]
    data$group = factor( as.character( data[ , colorBy ] ), levels = names( palette ) )

    # axis ticks in steps of 5
    xBreaks = seq( 5 * floor( min( data$UMAP_1 ) / 5 ), 5 * ceiling( max( data$UMAP_1 ) / 5 ), by = 5 )
    yBreaks = seq( 5 * floor( min( data$UMAP_2 ) / 5 ), 5 * ceiling( max( data$UMAP_2 ) / 5 ), by = 5 )
    xLimits = range( c( data$UMAP_1, xBreaks ) )
    yLimits = range( c( data$UMAP_2, yBreaks ) )

    p <- ggplot( data, aes( UMAP_1, UMAP_2, color = group ) ) + geom_point( size = diameterToSize( 1 ), alpha = 0.82, stroke = 0 ) +
      scale_color_manual( values = palette, drop = TRUE, name = NULL ) +
      scale_x_continuous( breaks = xBreaks, labels = cleanNumber, limits = xLimits, expand = expansion( mult = c( 0.025, 0.025 ) ) ) +
      scale_y_continuous( breaks = yBreaks, labels = cleanNumber, limits = yLimits, expand = expansion( mult = c( 0.025, 0.025 ) ) ) +
      coord_equal( clip = "off" ) + labs( x = "UMAP 1", y = "UMAP 2" ) +
      themePanel( box = TRUE ) +
      theme( legend.position = "inside", legend.position.inside = c( 0.5, 0.01 ), legend.justification = c( 0.5, 0 ), legend.direction = "horizontal",
             legend.text = element_text( size = labelSize, margin = margin( l = 1, r = 4 ) ), legend.key.size = unit( legendDotDiameter + 2, "pt" ),
             legend.key.spacing.x = unit( 0, "pt" ) ) +
      guides( color = guide_legend( override.aes = list( size = diameterToSize( legendDotDiameter ), alpha = 1 ), nrow = 1 ) )

    if( labelOnData )
    {
        # label at the median position of the cells of a cluster
        groups = unique( data$group )
        labelPos = data.frame( label = groups,
                               UMAP_1 = sapply( groups, function( g ) median( data$UMAP_1[ data$group == g ] ) ),
                               UMAP_2 = sapply( groups, function( g ) median( data$UMAP_2[ data$group == g ] ) ) )

        p <- p + guides( color = "none" ) +
          geom_text_repel( data = labelPos, aes( UMAP_1, UMAP_2, label = label ), inherit.aes = FALSE, family = "Arial", size = labelSize / .pt, color = "black",
                           min.segment.length = 0, segment.color = "black", segment.size = lw, box.padding = 0.2, point.padding = 0.16, max.overlaps = Inf, seed = 42 )
    }
    setPanelSize( p, umapWidth, umapWidth * diff( yLimits ) / diff( xLimits ) )
}


### Figure 3B: cell states (OPC and six oligodendrocyte clusters)
states = c( "OPC", paste0( "OL", 1:6 ) )
state_color = cluster_colors[ 1:length( states ) ]
names( state_color ) <- states

g <- plotUMAP( keptCells, "final_cluster_label", state_color, labelOnData = TRUE )
savePointsAsImage( g, paste0( figurePath, "Figure_3B_UMAP_OPC_OL_States" ) )

### Figure 3A: same cells colored by diagnosis
g <- plotUMAP( keptCells, "diagnosis", diagnosis_color )
savePointsAsImage( g, paste0( figurePath, "Figure_3A_UMAP_OPC_OL_Diagnosis" ) )

### all clusters before removing the small clusters and the clusters that come from few donors
cells$cluster_name <- paste0( "C", cells$natural_cluster )
allClusters = paste0( "C", sort( unique( cells$natural_cluster ) ) )
all_color = rep( cluster_colors, length.out = length( allClusters ) )
names( all_color ) <- allClusters

g <- plotUMAP( cells, "cluster_name", all_color, labelOnData = TRUE )
savePointsAsImage( g, paste0( figurePath, "UMAP_OPC_OL_AllClusters" ) )
