library(ggplot2)
library(gtable)
library(GO.db)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )

# output folder
outputPath = paste0( home_directory, "Figure_5/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )


units = c( "MSA_L1", "MSA_L3", "PD_L1", "PD_L3" )

# fitted curves of the genes that differ between disease and HC (made by 02_Disease_Dynamics_ConditionTest_Analysis.R)
# list with one entry per disease and lineage: genes, raw_response, log1p_response, scaled (z-score of log1p, genes x 100: 50 cell groups HC, 50 cell groups disease),
# time (median pseudotime of the cell groups, 0 = start and 1 = end of the compared range), median_pseudotime
curves <- readRDS( paste0( objectPath, "Disease_Dynamics_Curves.rds" ) )

# modules of the genes, overview of the modules and GO enrichment of the modules (made by 03_Disease_Dynamics_Modules_Analysis.R)
plotData <- readRDS( paste0( objectPath, "Disease_Dynamics_Modules_PlotData.rds" ) )
assignTab = plotData$assignTab
moduleTab = plotData$moduleTab
goTab = plotData$goTab


# range of the z-scores of a module (all curves of HC and disease), extended by 2% at both ends
moduleRange <- function( v )
{
    r = range( v )
    r + c( -1, 1 ) * 0.02 * diff( r )
}

hc_color = "#6FA8DC"
disease_color = "#E58F73"

# all curves as long table: one row per gene, condition and point; each point is drawn at the median pseudotime of its cell group
curveLong = c()
for( u in units )
{
    x = curves[[ u ]]
    genes = rownames( x$scaled )
    curveLong = rbind( curveLong, data.frame( unit_id = u, gene = rep( genes, each = 100 ), value = as.vector( t( x$scaled ) ),
                                              condition_segment = rep( rep( c( "HC", "Disease" ), each = 50 ), times = length( genes ) ),
                                              normalized_pseudotime = rep( rep( x$time, 2 ), times = length( genes ) ) ) )
}
cand = match( paste( curveLong$unit_id, curveLong$gene ), paste( assignTab$unit_id, assignTab$gene ) )
curveLong$disease = assignTab$disease[ cand ]
curveLong$assigned = assignTab$assigned[ cand ]
curveLong$module_id = assignTab$module_id[ cand ]
cand = match( curveLong$module_id, moduleTab$module_id )
curveLong$display_cluster = factor( moduleTab$display_cluster[ cand ], levels = paste0( "C", 1:max( table( moduleTab$unit_id ) ) ) )
curveLong$display_unit = factor( moduleTab$display_unit[ cand ], levels = gsub( "_", " ", units ) )
curveLong$condition = factor( ifelse( curveLong$condition_segment == "HC", "HC", curveLong$disease ), levels = c( "HC", "MSA", "PD" ) )
curveLong$display_condition = factor( curveLong$condition_segment, levels = c( "HC", "Disease" ) )
curveLong = curveLong[ curveLong$assigned, ]

# mean curve of each module
mainCurves = unique( curveLong[ , c( "unit_id", "module_id", "display_cluster", "display_unit", "condition", "display_condition", "normalized_pseudotime" ) ] )
curveKey = paste( curveLong$module_id, curveLong$display_condition, curveLong$normalized_pseudotime )
meanValue = tapply( curveLong$value, curveKey, mean )
mainCurves$main = as.numeric( meanValue[ paste( mainCurves$module_id, mainCurves$display_condition, mainCurves$normalized_pseudotime ) ] )

# Figure 5A shows each module over its own range (no axis): z-scores scaled to 0 - 1 within the range of the module
ranges = lapply( split( curveLong$value, curveLong$module_id ), moduleRange )
lower = sapply( ranges, "[", 1 )
width = sapply( ranges, diff )
curveLong$value_panel = ( curveLong$value - lower[ curveLong$module_id ] ) / width[ curveLong$module_id ]
mainCurves$main_panel = ( mainCurves$main - lower[ mainCurves$module_id ] ) / width[ mainCurves$module_id ]

panelTab = data.frame( moduleTab[ , c( "unit_id", "module_id", "genes" ) ],
                       display_cluster = factor( moduleTab$display_cluster, levels = levels( curveLong$display_cluster ) ),
                       display_unit = factor( moduleTab$display_unit, levels = levels( curveLong$display_unit ) ) )

### Figure 5A: all modules of both diseases and both lineages
# size of one cell of the grid in points
overviewCell = 49.3

p_all <- ggplot() +
  geom_line( data = curveLong, aes( normalized_pseudotime, value_panel, colour = display_condition, group = interaction( unit_id, module_id, gene, display_condition ) ), linewidth = lw, alpha = 0.18, lineend = "round" ) +
  geom_line( data = mainCurves, aes( normalized_pseudotime, main_panel, colour = display_condition, group = interaction( unit_id, module_id, display_condition ) ), linewidth = lw, lineend = "round" ) +
  geom_text( data = panelTab, aes( x = 0.96, y = 0.96, label = paste0( "n = ", genes ) ), inherit.aes = FALSE, hjust = 1, vjust = 1, family = "Arial",
             size = annotationSize / .pt, colour = "black" ) +
  facet_grid( display_unit ~ display_cluster, switch = "y", drop = FALSE ) +
  scale_colour_manual( values = c( HC = hc_color, Disease = disease_color ), breaks = c( "HC", "Disease" ), drop = FALSE ) +
  scale_x_continuous( expand = c( 0, 0 ) ) + scale_y_continuous( expand = c( 0, 0 ) ) + coord_fixed( ratio = 1, xlim = c( 0, 1 ), ylim = c( 0, 1 ), expand = FALSE, clip = "on" ) +
  guides( colour = guide_legend( nrow = 1, override.aes = list( linewidth = lw, alpha = 1 ) ) ) +
  theme_classic( base_family = "Arial", base_size = titleSize ) +
  theme( text = element_text( family = "Arial", size = titleSize, colour = "black" ), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
         panel.border = element_rect( fill = NA, colour = "black", linewidth = lw ), panel.spacing.x = unit( 0, "pt" ), panel.spacing.y = unit( 0, "pt" ),
         strip.placement = "outside", strip.background.x = element_blank(), strip.text.x = element_text( family = "Arial", size = titleSize, margin = margin( 0, 0, 2, 0 ) ),
         strip.background.y = element_rect( fill = stripFill, colour = "black", linewidth = lw ),
         strip.text.y.left = element_text( family = "Arial", angle = 90, size = titleSize, margin = margin( 2, 3, 2, 3 ) ),
         legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text( family = "Arial", size = labelSize ),
         legend.key.width = unit( 18, "pt" ), legend.key = element_blank(), plot.margin = margin( 4, 3, 2, 4 ) )

saveGrob( setPanelSize( p_all, overviewCell, overviewCell ), paste0( figurePath, "Figure_5A_Modules_Overview" ) )


#
# Plot of one module: curves of all genes with the mean curve (left) and the most significant GO terms (right)
#

# size of the plotting areas in points: curves (width = height), GO terms (at least this width, the height of the curves)
curveSize = 80
goWidth = 150
goTextPadding = 3

# names of the GO terms: official name in GO.db, first letter upper case
goSets = unique( as.data.frame( readRDS( paste0( dataPath, "ReferenceFiles/GeneSets/MSigDB_2026.1_GO_BP.rds" ) ) )[ , c( "gs_id", "gs_exact_source" ) ] )
goTermName <- function( termId )
{
    x = unname( Term( goSets$gs_exact_source[ match( termId, goSets$gs_id ) ] ) )
    paste0( toupper( substr( x, 1, 1 ) ), substring( x, 2 ) )
}
textWidth <- function( x, size ) withPdfDevice( max( sapply( x, function( v ) convertWidth( grobWidth( textGrob( v, gp = gpar( fontfamily = "Arial", fontsize = size ) ) ), "pt", valueOnly = TRUE ) ) ) )

moduleTheme <- function( showAxes = TRUE, showLegend = TRUE )
{
    theme_classic( base_family = "Arial", base_size = titleSize ) +
      theme( text = element_text( family = "Arial", size = titleSize, colour = "black" ),
             plot.title = element_text( family = "Arial", size = titleSize, hjust = 0.5, margin = margin( 0, 0, 3, 0 ) ),
             axis.title = if( showAxes ) element_text( family = "Arial", size = axisTitleSize ) else element_blank(),
             axis.text = if( showAxes ) element_text( family = "Arial", size = tickTextSize, colour = "black" ) else element_blank(),
             axis.ticks = if( showAxes ) element_line( linewidth = lw, colour = "black" ) else element_blank(),
             axis.ticks.length = unit( tickLength, "pt" ), axis.line = element_blank(), panel.border = element_rect( fill = NA, colour = "black", linewidth = lw ), panel.grid = element_blank(),
             legend.position = if( showLegend ) "bottom" else "none", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text( family = "Arial", size = labelSize ),
             legend.key.width = unit( 18, "pt" ), legend.key = element_blank(), legend.margin = margin( 0, 0, 0, 0 ), plot.margin = margin( 3, 3, 3, 3 ) ) +
      theme( strip.background = element_blank(), strip.text = element_text( family = "Arial", size = titleSize ), panel.spacing = unit( 4, "pt" ) )
}

fdrStars <- function( x ) ifelse( x < 0.001, "***", ifelse( x < 0.01, "**", ifelse( x < 0.05, "*", "" ) ) )

plotModule <- function( mid )
{
    meta = moduleTab[ moduleTab$module_id == mid, ]
    sub = curveLong[ curveLong$module_id == mid, ]
    main = mainCurves[ mainCurves$module_id == mid, ]
    diseaseLabel = unique( as.character( sub$disease ) )
    statPt = annotationSize / .pt

    # y axis over the range of the module, ticks at even numbers (every 4 if the range is wider than 8)
    yLim = moduleRange( sub$value )
    yTop = yLim[ 2 ] - 0.04 * diff( yLim )
    yStep = if( diff( yLim ) > 8 ) 4 else 2
    yBreaks = seq( ceiling( yLim[ 1 ] / yStep ) * yStep, floor( yLim[ 2 ] / yStep ) * yStep, by = yStep )

    p_curve <- ggplot() +
      geom_line( data = sub, aes( normalized_pseudotime, value, colour = condition, group = interaction( gene, condition ) ), linewidth = lw, alpha = 0.18, lineend = "round" ) +
      geom_line( data = main, aes( normalized_pseudotime, main, colour = condition, group = condition ), linewidth = lw, lineend = "round" )
    # modules that are enriched in the HOG ranking get the NES written into the plot, the number of genes is below it
    if( isTRUE( meta$HOG_FDR05 ) )
    {
        p_curve <- p_curve + annotate( "text", x = 0.04, y = yTop, label = sprintf( "NES = %.2f%s", meta$HOG_NES, fdrStars( meta$HOG_FDR ) ), hjust = 0, vjust = 1, family = "Arial",
                                       size = statPt, colour = "black" ) +
          annotate( "text", x = 0.04, y = yTop - 0.12 * diff( yLim ), label = paste0( meta$genes, " genes" ), hjust = 0, vjust = 1, family = "Arial", size = statPt, colour = "black" )
    }
    else
    {
        p_curve <- p_curve + annotate( "text", x = 0.96, y = yTop, label = paste0( meta$genes, " genes" ), hjust = 1, vjust = 1, family = "Arial", size = statPt, colour = "black" )
    }
    p_curve <- p_curve +
      scale_colour_manual( values = setNames( c( hc_color, disease_color ), c( "HC", diseaseLabel ) ), breaks = c( "HC", diseaseLabel ), drop = FALSE ) +
      scale_x_continuous( breaks = c( 0, 0.5, 1 ), labels = c( "0", "0.5", "1" ), expand = c( 0, 0 ) ) +
      scale_y_continuous( breaks = yBreaks, labels = as.character( yBreaks ), expand = c( 0, 0 ) ) + coord_fixed( ratio = 1 / diff( yLim ), xlim = c( 0, 1 ), ylim = yLim, expand = FALSE, clip = "on" ) +
      labs( title = paste( gsub( "_", " ", meta$unit_id ), meta$display_cluster ), x = "OL pseudotime", y = "Z-score" ) +
      guides( colour = guide_legend( nrow = 1, override.aes = list( linewidth = lw, alpha = 1 ) ) ) + moduleTheme( showAxes = TRUE, showLegend = TRUE )

    #the 6 most significant GO terms as bars
    go = goTab[ goTab$module_id == mid & goTab$FDR05, ]
    go = go[ order( go$padj_within_module_BH, go$pvalue ), ]
    width = goWidth
    if( nrow( go ) == 0 )
    {
        p_go <- ggplot() + annotate( "text", x = 0.5, y = 0.5, label = "No significant GO BP", family = "Arial", size = labelSize / .pt, colour = "#6F6F6F" ) +
          scale_x_continuous( limits = c( 0, 1 ), expand = c( 0, 0 ) ) + scale_y_continuous( limits = c( 0, 1 ), expand = c( 0, 0 ) ) + labs( title = "GO BP" ) +
          moduleTheme( showAxes = FALSE, showLegend = FALSE )
    }
    else
    {
        go = go[ 1:min( 6, nrow( go ) ), ]
        go$label = goTermName( go$term_id )
        # the plotting area is wider if the longest name does not fit into it
        width = max( goWidth, textWidth( go$label, barTextSize ) / 0.982 + goTextPadding )

        go$score = -log10( pmax( go$padj_within_module_BH, .Machine$double.xmin ) )
        go$row = factor( 1:nrow( go ), levels = rev( 1:nrow( go ) ) )
        xmax = max( go$score ) * 1.03

        tickNumber <- function( x )
        {
            res = trimws( formatC( x, format = "fg", digits = 4 ) )
            res = sub( "([.][0-9]*[1-9])0+$", "\\1", res )
            res = sub( "[.]0+$", "", res )
            res[ x == 0 ] = "0"
            res
        }

        p_go <- ggplot( go, aes( score, row ) ) + geom_col( fill = "#D9D9D9", width = 0.86 ) +
          geom_text( aes( x = xmax * 0.018, label = label ), hjust = 0, family = "Arial", size = barTextSize / .pt, colour = "black" ) +
          scale_x_continuous( limits = c( 0, xmax ), labels = tickNumber, expand = c( 0, 0 ) ) + scale_y_discrete( labels = NULL ) +
          labs( title = "GO BP", x = expression( paste( "−log"[10], " FDR" ) ), y = NULL ) + moduleTheme( showAxes = TRUE, showLegend = FALSE ) +
          theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin( 3, 8, 3, 3 ) )
    }

    #curves and GO terms next to each other
    g_curve = setPanelSize( p_curve, curveSize, curveSize )
    g_go = setPanelSize( p_go, width, curveSize )
    # both panels at the same height: the parts get empty space above and below up to the larger title and the larger part below the panel (legend of the curves)
    panelTop <- function( g ) withPdfDevice( convertHeight( sum( g$heights[ seq_len( min( g$layout$t[ grepl( "^panel", g$layout$name ) ] ) - 1 ) ] ), "pt", valueOnly = TRUE ) )
    tableHeight <- function( g ) withPdfDevice( convertHeight( sum( g$heights ), "pt", valueOnly = TRUE ) )
    top = max( panelTop( g_curve ), panelTop( g_go ) )
    bottom = max( tableHeight( g_curve ) - panelTop( g_curve ), tableHeight( g_go ) - panelTop( g_go ) )
    padRows <- function( g )
    {
        g = gtable_add_rows( g, unit( top - panelTop( g ), "pt" ), pos = 0 )
        gtable_add_rows( g, unit( top + bottom - tableHeight( g ), "pt" ), pos = -1 )
    }
    g_curve = padRows( g_curve )
    g_go = padRows( g_go )
    fig = gtable( widths = unit.c( sum( g_curve$widths ), unit( 4, "pt" ), sum( g_go$widths ) ), heights = sum( g_curve$heights ) )
    gtable_add_grob( fig, list( g_curve, g_go ), t = 1, l = c( 1, 3 ) )
}

### one figure per module (Figure 5B, C, F and G are PD L1 C1, PD L3 C1, MSA L1 C3 and MSA L3 C1)
for( k in 1:nrow( moduleTab ) )
{
    modulePath = paste0( figurePath, "Modules_GO_BP/", moduleTab$unit_id[ k ], "/" )
    dir.create( modulePath, recursive = TRUE, showWarnings = FALSE )
    saveGrob( plotModule( moduleTab$module_id[ k ] ), paste0( modulePath, moduleTab$unit_id[ k ], "_", moduleTab$display_cluster[ k ] ) )
}

# which figure file is which panel of Figure 5 (file names relative to the figure folder), 06_SelectedGenes_OXPHOS_PD_Plots.R adds panel D and E
panelFiles = data.frame( panel = c( "A", "B", "C", "F", "G" ),
                         disease = c( "MSA and PD", "PD", "PD", "MSA", "MSA" ),
                         lineage = c( "all", "L1", "L3", "L1", "L3" ),
                         module_id = c( "all modules", "PD_L1__D01", "PD_L3__D01", "MSA_L1__D03", "MSA_L3__D01" ),
                         figure_pdf = c( "Figure_5A_Modules_Overview.pdf", "Modules_GO_BP/PD_L1/PD_L1_C1.pdf", "Modules_GO_BP/PD_L3/PD_L3_C1.pdf", "Modules_GO_BP/MSA_L1/MSA_L1_C3.pdf", "Modules_GO_BP/MSA_L3/MSA_L3_C1.pdf" ) )
write.table( panelFiles, file = paste0( figurePath, "Figure_5_PanelFiles.txt" ), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t" )
