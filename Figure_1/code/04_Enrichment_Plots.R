library(ggplot2)
library(stringr)
library(scales)
library(gtable)
library(GO.db)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_1/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )
tablePath = paste0( outputPath, "tables/" )
objectPath = paste0( outputPath, "objects/" )

source( paste0( codePath, "saveFigure.R" ) )
source( paste0( codePath, "panelLayout.R" ) )

writeTable <- function( x, file, quote = FALSE ) write.table( x, file = paste0( tablePath, file ), col.names = TRUE, row.names = FALSE, quote = quote, sep = "\t" )


# size of the plotting areas in points and titles of the panels
panelSize <- list( KEGG = c( 63.4, 192.2 ), Reactome = c( 58.6, 158.6 ), GO_BP = c( 61.3, 186 ), GO_CC = c( 58.8, 72.4 ), GO_MF = c( 49.7, 152.1 ) )
panelTitle <- c( KEGG = "KEGG Enrichment", Reactome = "Reactome Enrichment", GO_BP = "GO-BP Enrichment", GO_CC = "GO-CC Enrichment", GO_MF = "GO-MF Enrichment" )

# line width of the circles in the unit of the stroke argument of geom_point
strokeFromPoint <- function( pt ) pt * 96 / 72 * 2 / .stroke

# names of the gene sets: official name (GO terms: name in GO.db), first letter upper case, two words per line
goIds <- unlist( lapply( c( "GO_BP", "GO_CC", "GO_MF" ), function( db ) {
    sets = unique( as.data.frame( readRDS( paste0( dataPath, "ReferenceFiles/GeneSets/MSigDB_2026.1_", db, ".rds" ) ) )[ , c( "gs_name", "gs_exact_source" ) ] )
    setNames( sets$gs_exact_source, sets$gs_name ) } ) )
formatTermName <- function( id, description )
{
    x = ifelse( id %in% names( goIds ), unname( Term( goIds[ id ] ) ), description )
    x = paste0( toupper( substr( x, 1, 1 ) ), substring( x, 2 ) )
    vapply( strsplit( x, "\\s+" ), function( w ) paste( sapply( split( w, ceiling( seq_along( w ) / 2 ) ), paste, collapse = " " ), collapse = "\n" ), "" )
}

# the names of the gene sets get new line breaks so that every name needs at most two lines
rewrapTerms <- function( terms )
{
    flat = gsub( "\\s+", " ", terms )
    for( w in 14:60 )
    {
        res = vapply( flat, function( s ) paste( strwrap( s, width = w ), collapse = "\n" ), "" )
        if( max( lengths( strsplit( res, "\n" ) ) ) <= 2 ) return( res )
    }
    res
}

plotEnrichment <- function( plotTab, outFile, db )
{
    plotTab$display_term = rewrapTerms( plotTab$display_term )
    cmin = min( plotTab$Count )
    cmax = max( plotTab$Count )
    sizeBreaks = if( cmin == cmax ) cmin else unique( round( seq( cmin, cmax, length.out = 3 ) ) )

    fdrLabel <- function( x ) sapply( x, function( v ) if( !is.finite( v ) ) NA else if( abs( v ) < .Machine$double.eps ) "0" else scientific_format( digits = 1 )( v ) )

    p <- ggplot( plotTab, aes( x = gene_ratio, y = y_pos ) ) +
      geom_point( aes( size = Count, fill = p.adjust ), shape = 21, color = "grey30", stroke = strokeFromPoint( lineWidth ) ) +
      scale_fill_gradient( low = "#E06666", high = "#6FA8DC", limits = c( 0, 0.05 ), breaks = c( 0, 0.0125, 0.025, 0.0375, 0.05 ),
                           labels = fdrLabel, name = "Padj",
                           guide = guide_colorbar( order = 2, barwidth = unit( 6, "pt" ), barheight = unit( 50, "pt" ), frame.colour = "black", frame.linewidth = lw, ticks.colour = "black", ticks.linewidth = lw ) ) +
      scale_size_continuous( name = "Count", range = diameterToSize( c( 4.8, 12.5 ) ), breaks = sizeBreaks, limits = c( cmin, cmax ),
                             guide = guide_legend( order = 1, override.aes = list( fill = "grey50" ) ) ) +
      scale_x_continuous( expand = expansion( mult = c( 0.12, 0.2 ) ), breaks = function( lim ) seq( lim[ 1 ], lim[ 2 ], length.out = 3 ),
                          labels = number_format( accuracy = 0.01 ) ) +
      scale_y_continuous( breaks = plotTab$y_pos, labels = plotTab$display_term,
                          expand = expansion( add = c( max( plotTab$row_height ) / 2 + 0.05, max( plotTab$row_height ) / 2 + 0.05 ) ) ) +
      coord_cartesian( clip = "off" ) + labs( x = "Gene Ratio", y = NULL, title = panelTitle[ db ] ) +
      themePanel( box = TRUE ) +
      theme( axis.text.y = element_text( size = labelSize, lineheight = 0.8, color = "black" ), axis.title.y = element_blank(),
             plot.title.position = "panel",
             panel.grid.major = element_line( color = "#F0F0F0", linewidth = lw ), panel.grid.minor = element_blank(),
             legend.position = "right", legend.justification = "top", legend.key.size = unit( 14, "pt" ), legend.spacing.y = unit( 6, "pt" ),
             legend.title = element_text( size = axisTitleSize, margin = margin( b = 3 ) ), legend.text = element_text( size = legendValueSize ) )

    saveGrob( setPanelSize( p, panelSize[[ db ]][ 1 ], panelSize[[ db ]][ 2 ] ), outFile )
}


# plot of the 8 most significant gene sets of every collection, enrichment of all differentially expressed genes together (03_Enrichment_Analysis.R)
resultTables = readRDS( paste0( objectPath, "Enrichment_Results.rds" ) )
for( db in c( "GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome" ) )
{
    plotTab = resultTables[[ db ]][ resultTables[[ db ]]$p.adjust < 0.05, ]
    if( nrow( plotTab ) == 0 ) next

    ratio = str_split_fixed( plotTab$GeneRatio, "/", 2 )
    plotTab$gene_ratio = as.numeric( ratio[ , 1 ] ) / as.numeric( ratio[ , 2 ] )
    plotTab = plotTab[ order( plotTab$p.adjust, -plotTab$Count, -plotTab$gene_ratio ), ]
    plotTab = plotTab[ 1:min( 8, nrow( plotTab ) ), ]
    plotTab$display_term = formatTermName( plotTab$ID, plotTab$Description )
    plotTab$neg_log10_fdr = -log10( plotTab$p.adjust )
    # most significant gene set on top
    plotTab = plotTab[ order( plotTab$neg_log10_fdr ), ]
    plotTab$row_height = 3.28
    plotTab$y_pos = cumsum( plotTab$row_height ) - plotTab$row_height / 2

    # quotes are needed here because of the line breaks in the names
    writeTable( plotTab, paste0( "Enrichment/Enrichment_", db, "_All_DEG_PlottedTerms.txt" ), quote = TRUE )
    plotEnrichment( plotTab, paste0( figurePath, "Figure_1_Enrichment_", db ), db )
}
