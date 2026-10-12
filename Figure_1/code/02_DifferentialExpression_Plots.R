library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(scales)

# change to your project directory
home_directory <- "path/to/project_directory/"

outputPath = paste0( home_directory, "Figure_1/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

# results of 01_DifferentialExpression_Analysis.R
dds <- readRDS( paste0( objectPath, "HOG_DESeq2_dds.rds" ) )
res <- readRDS( paste0( objectPath, "HOG_DESeq2_results.rds" ) )
rld <- readRDS( paste0( objectPath, "HOG_rlog.rds" ) )
plotData <- readRDS( paste0( objectPath, "HOG_PlotData.rds" ) )
df_pca = plotData$pca

# axis labels: write 0 and not 0.0 or 0.00
# numbers without a zero at the end when they have decimals (0.50 is written 0.5)
zeroLabel <- function( x )
{
    lab = label_number()( x )
    lab = ifelse( grepl( "\\.", lab ), sub( "\\.?0+$", "", lab ), lab )
    ifelse( x == 0, "0", lab )
}

condition_color = c( Control = "#6FA8DC", asyn_PFF = "#E06666" )

qc_theme <- theme_classic( base_family = "Arial", base_size = 12 ) +
  theme( text = element_text( family = "Arial", face = "bold" ),
         axis.text = element_text( face = "bold" ), axis.title = element_text( face = "bold" ),
         plot.title = element_text( face = "bold", size = 13 ),
         panel.grid.major = element_line( color = "#F0F0F0" ), panel.grid.minor = element_blank(),
         panel.border = element_rect( color = "black", fill = NA, linewidth = 0.5 ),
         legend.position = "right" )

#
### Quality control
#

### library size in million reads
libSize <- colSums( counts( dds ) )
df_lib = data.frame( sample = names( libSize ), libsize = libSize / 1e6, cond = dds$condition )
p <- ggplot( df_lib, aes( sample, libsize, fill = cond ) ) +
  geom_col( width = 0.75, color = "black", linewidth = 0.3 ) +
  scale_fill_manual( values = condition_color, name = "condition" ) +
  labs( x = NULL, y = "Library size (million reads)" ) + qc_theme +
  theme( axis.text.x = element_text( angle = 45, hjust = 1, face = "bold" ), panel.grid.major.x = element_blank() ) +
  scale_y_continuous( labels = zeroLabel )
saveFigure( p, paste0( figurePath, "QC/QC_LibrarySize" ), width = 6, height = 4 )

### PCA
p <- ggplot( df_pca, aes( PC1, PC2, fill = cond, label = sample ) ) +
  geom_point( size = 4, shape = 21, color = "black", stroke = 0.5 ) +
  geom_text_repel( size = 3, family = "Arial", fontface = "bold", max.overlaps = Inf, seed = 42 ) +
  scale_fill_manual( values = condition_color, name = "condition" ) +
  labs( x = sprintf( "PC1 (%.1f%%)", df_pca$PC1_var_pct[ 1 ] ), y = sprintf( "PC2 (%.1f%%)", df_pca$PC2_var_pct[ 1 ] ) ) + qc_theme +
  scale_x_continuous( labels = zeroLabel ) + scale_y_continuous( labels = zeroLabel )
saveFigure( p, paste0( figurePath, "QC/QC_PCA" ), width = 5.5, height = 4.5 )

### heatmap of the distances between samples
sampleDist <- as.matrix( dist( t( assay( rld ) ) ) )
upper <- ceiling( max( sampleDist ) / 10 ) * 10
dist_color <- colorRamp2( c( 0, upper / 2, upper ), c( "#FFFFFF", "#F4B6B6", "#E06666" ) )
ha <- HeatmapAnnotation( Condition = as.character( dds$condition ),
                         col = list( Condition = condition_color ),
                         annotation_name_gp = gpar( fontfamily = "Arial", fontface = "bold", fontsize = 10 ) )
hm <- Heatmap( sampleDist, col = dist_color, name = "Euclidean", top_annotation = ha,
               width = unit( 52, "mm" ), height = unit( 52, "mm" ),
               rect_gp = gpar( col = "#E0E0E0", lwd = 0.4 ), border_gp = gpar( col = "black", lwd = 0.8 ),
               row_names_gp = gpar( fontfamily = "Arial", fontface = "bold", fontsize = 10 ),
               column_names_gp = gpar( fontfamily = "Arial", fontface = "bold", fontsize = 10 ),
               heatmap_legend_param = list( at = seq( 0, upper, length.out = 5 ),
                                            labels = format( seq( 0, upper, length.out = 5 ), trim = TRUE ),
                                            title_gp = gpar( fontfamily = "Arial", fontface = "bold", fontsize = 10 ),
                                            labels_gp = gpar( fontfamily = "Arial", fontface = "bold", fontsize = 9 ) ) )
saveFigure( function() draw( hm ), paste0( figurePath, "QC/QC_SampleDistanceHeatmap" ), width = 5, height = 5 )

### dispersion estimates
saveFigure( function() {
  par( family = "Arial", font = 2, font.axis = 2, font.lab = 2, font.main = 2 )
  plotDispEsts( dds, main = NULL ) }, paste0( figurePath, "QC/QC_DispersionEstimates" ), width = 5.5, height = 4.5 )

### histogram of p-values
df_res <- as.data.frame( res )
p <- ggplot( df_res[ !is.na( df_res$padj ) & !is.na( df_res$pvalue ), ], aes( pvalue ) ) +
  geom_histogram( bins = 40, fill = "#6FA8DC", color = "black", linewidth = 0.3 ) +
  geom_vline( xintercept = 0.05, linetype = "dashed", color = "#E06666", linewidth = 0.5 ) +
  labs( x = "raw p-value", y = "genes" ) + qc_theme + theme( panel.grid.major.x = element_blank() ) +
  scale_x_continuous( labels = zeroLabel ) + scale_y_continuous( labels = zeroLabel )
saveFigure( p, paste0( figurePath, "QC/QC_PValueHistogram" ), width = 5, height = 4 )

### MA plot
df_ma <- df_res[ !is.na( df_res$padj ), ]
df_ma$sig <- factor( ifelse( df_ma$padj < 0.05, ifelse( df_ma$log2FoldChange > 0, "Up", "Down" ), "NS" ), levels = c( "Up", "Down", "NS" ) )
p <- ggplot( df_ma, aes( log2( baseMean + 1 ), log2FoldChange, color = sig ) ) +
  geom_point( alpha = 0.55, size = 0.6 ) + geom_hline( yintercept = 0, color = "black", linewidth = 0.3 ) +
  scale_color_manual( values = c( Up = "#E06666", Down = "#6FA8DC", NS = "#E0E0E0" ), breaks = c( "Up", "Down", "NS" ), name = "padj < 0.05" ) +
  labs( x = expression( log[2] * "(baseMean + 1)" ), y = expression( log[2] * "(FoldChange (MLE))" ) ) +
  qc_theme + guides( color = guide_legend( override.aes = list( size = 2, alpha = 1 ) ) ) +
  scale_x_continuous( labels = zeroLabel ) + scale_y_continuous( labels = zeroLabel )
saveFigure( p, paste0( figurePath, "QC/QC_MAPlot" ), width = 6, height = 4.5 )

#
### Volcano plot (Figure 1C)
#

volcanoData <- plotData$volcano
labelGenes <- plotData$labels
n_up <- sum( volcanoData$significance == "Up" )
n_down <- sum( volcanoData$significance == "Down" )

# x-axis: symmetric, ticks at round numbers
xMax <- max( abs( volcanoData$log2FoldChange ) )
xPad <- max( 0.35, xMax * 0.08 )
xRange <- c( -xMax - xPad, xMax + xPad )
tickStep <- max( 0.5, signif( xMax / 2, 1 ) )
xTicks <- seq( -floor( xMax / tickStep ) * tickStep, floor( xMax / tickStep ) * tickStep, by = tickStep )
yMax <- max( volcanoData$logFDR ) * 1.08

# sizes of the volcano plot, in points
volcanoPointDiameter = 2.4
volcanoPanelWidth = 155.9
volcanoPanelHeight = 146.4

p <- ggplot( volcanoData, aes( x = log2FoldChange, y = logFDR, color = significance ) ) +
  geom_point( alpha = 0.62, size = diameterToSize( volcanoPointDiameter ), stroke = 0 ) +
  geom_hline( yintercept = -log10( 0.05 ), linetype = "dashed", color = "#7F8C8D", linewidth = lw ) +
  scale_color_manual( values = c( Up = "#E06666", Down = "#6FA8DC", NS = "#E0E0E0" ),
                      labels = c( Up = sprintf( "Up (%d)", n_up ), Down = sprintf( "Down (%d)", n_down ), NS = "NS" ),
                      breaks = c( "Up", "Down", "NS" ) ) +
  annotate( "text", x = xRange[ 2 ] - diff( xRange ) * 0.035, y = yMax * 0.94, label = sprintf( "Up %d", n_up ),
            hjust = 1, vjust = 1, family = "Arial", size = annotationSize / .pt, color = "#FF0000" ) +
  annotate( "text", x = xRange[ 1 ] + diff( xRange ) * 0.035, y = yMax * 0.94, label = sprintf( "Down %d", n_down ),
            hjust = 0, vjust = 1, family = "Arial", size = annotationSize / .pt, color = "#0000FF" ) +
  geom_text_repel( data = labelGenes, aes( label = symbol ), size = labelSize / .pt, max.overlaps = 20,
                   box.padding = 0.25, point.padding = 0.12, segment.color = "grey40", segment.size = lw,
                   color = "black", fontface = "italic", min.segment.length = 0, family = "Arial",
                   xlim = xRange, ylim = c( 0, yMax ) ) +
  coord_cartesian( xlim = xRange, ylim = c( 0, yMax ), expand = FALSE ) +
  scale_x_continuous( breaks = xTicks, labels = zeroLabel ) + scale_y_continuous( labels = zeroLabel ) +
  labs( x = expression( log[2] * FC ), y = expression( paste( "−log"[10], " FDR" ) ) ) +
  themePanel( box = TRUE ) +
  theme( panel.grid.major = element_line( color = "#F0F0F0", linewidth = lw ), panel.grid.minor = element_blank(),
         legend.position = "none" )
saveGrob( setPanelSize( p, volcanoPanelWidth, volcanoPanelHeight ), paste0( figurePath, "Figure_1C_VolcanoPlot" ) )
