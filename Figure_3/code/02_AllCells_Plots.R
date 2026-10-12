library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_3/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )


# cell annotation of the snRNA-seq data set (columns: cell_id, sample_id, diagnosis, labels_major, labels_minor, UMAP_1, UMAP_2)
# made by DataPreparation/snRNAseq_HC_MSA_PD/code/CellAnnotation_AuthorLabels_UMAP.py from the data of Rydbirk et al. (2026)
cells = as.data.frame( fread( paste0( dataPath, "snRNAseq_HC_MSA_PD/NatComm_CellAnnotation.txt.gz" ) ) )

# cells are drawn in the order of the table, sort them by cell type
# (method radix: same order on every computer, independent of the language settings)
cells <- cells[ order( cells$labels_major, cells$labels_minor, cells$diagnosis, cells$sample_id, cells$cell_id, method = "radix" ), ]
cells$diagnosis <- factor( cells$diagnosis, levels = c( "HC", "MSA", "PD" ) )

# colors
diagnosis_color = c( HC = "#D0CDD6", MSA = "#E06666", PD = "#6FA8DC" )
major_color = c( Oligodendrocytes = "#E06666", OPCs = "#F2CC5E", Astrocytes = "#95D5B2", Microglia = "#B39BE8", `Medium spiny neurons` = "#3B86C4",
                 `GABAergic neurons` = "#7A5CCF", `GLUergic neurons` = "#4AA6A6", `Pericytes/endothelial` = "#5F6B77", PVMs = "#D9826B", Immune = "#5DA271" )
minor_color = c( OL_LINC01608 = "#D75D76", OL_SLC5A11 = "#E78DA0", OL_SGCZ = "#A94E76", OPCs = "#F2CC5E", COPs = "#C9A75D",
                 AS_homeostatic = "#95D5B2", AS_reactive = "#5DA271", MSN = "#3B86C4", `Inh. neurons` = "#7A5CCF", `Exc. neurons` = "#4AA6A6",
                 `MIC_steady-state` = "#C7B7EC", MIC_intermediate1 = "#A88EDB", MIC_intermediate2 = "#8766C8", MIC_activated = "#694BB0",
                 `Pericytes/endothelial` = "#5F6B77", PVMs = "#D9826B", Immune = "#5DA271" )


# width of the plotting area in points, the height follows from the same scale on both axes
umapWidth = 200
legendDotDiameter = 6


# UMAP colored by one column of the cell annotation
# labelOnData: write the names into the plot instead of a legend; a legend of up to three groups is in one row at the bottom of the plotting area, a longer one right of it
plotUMAP <- function( data, colorBy, palette, labelOnData = FALSE, pointSize = diameterToSize( 1 ) )
{
    palette = palette[ names( palette ) %in% unique( as.character( data[ , colorBy ] ) ) ]
    data$group = factor( as.character( data[ , colorBy ] ), levels = names( palette ) )

    # axis ticks in steps of 5
    xBreaks = seq( 5 * floor( min( data$UMAP_1 ) / 5 ), 5 * ceiling( max( data$UMAP_1 ) / 5 ), by = 5 )
    yBreaks = seq( 5 * floor( min( data$UMAP_2 ) / 5 ), 5 * ceiling( max( data$UMAP_2 ) / 5 ), by = 5 )
    xLimits = range( c( data$UMAP_1, xBreaks ) )
    yLimits = range( c( data$UMAP_2, yBreaks ) )

    p <- ggplot( data, aes( UMAP_1, UMAP_2, color = group ) ) + geom_point( size = pointSize, alpha = 0.82, stroke = 0 ) +
      scale_color_manual( values = palette, drop = TRUE, name = NULL ) +
      scale_x_continuous( breaks = xBreaks, labels = cleanNumber, limits = xLimits, expand = expansion( mult = c( 0.025, 0.025 ) ) ) +
      scale_y_continuous( breaks = yBreaks, labels = cleanNumber, limits = yLimits, expand = expansion( mult = c( 0.025, 0.025 ) ) ) +
      coord_equal( clip = "off" ) + labs( x = "UMAP 1", y = "UMAP 2" ) +
      themePanel( box = TRUE ) +
      theme( legend.text = element_text( size = labelSize, margin = margin( l = 1, r = 4 ) ), legend.key.size = unit( legendDotDiameter + 2, "pt" ),
             legend.key.spacing.x = unit( 0, "pt" ), legend.key.spacing.y = unit( 0, "pt" ) )
    if( length( palette ) <= 3 )
    {
        p <- p + theme( legend.position = "inside", legend.position.inside = c( 0.5, 0.01 ), legend.justification = c( 0.5, 0 ), legend.direction = "horizontal" ) +
          guides( color = guide_legend( override.aes = list( size = diameterToSize( legendDotDiameter ), alpha = 1 ), nrow = 1 ) )
    }
    else
    {
        p <- p + theme( legend.position = "right", legend.direction = "vertical" ) +
          guides( color = guide_legend( override.aes = list( size = diameterToSize( legendDotDiameter ), alpha = 1 ), ncol = 1 ) )
    }

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



### major cell types
g <- plotUMAP( cells, "labels_major", major_color, labelOnData = TRUE )
savePointsAsImage( g, paste0( figurePath, "UMAP_AllCells_CellTypes" ) )

### subtypes of the authors
g <- plotUMAP( cells, "labels_minor", minor_color )
savePointsAsImage( g, paste0( figurePath, "UMAP_AllCells_CellSubtypes" ) )

### diagnosis
g <- plotUMAP( cells, "diagnosis", diagnosis_color )
savePointsAsImage( g, paste0( figurePath, "UMAP_AllCells_Diagnosis" ) )
