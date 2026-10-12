library(clusterProfiler)
library(SingleCellExperiment)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/" )

# output folder
outputPath = paste0( home_directory, "Figure_4/" )
tablePath = paste0( outputPath, "tables/HC_Dynamics/" )
objectPath = paste0( outputPath, "objects/" )


lineageTab = data.frame( unit_id = c( "HC_L1_TDE", "HC_L2_TDE", "HC_L3_TDE" ), trajectory_id = c( "L1", "L2", "L3" ), model_lineage = 1:3, endpoint = c( "OL5", "OL6", "OL3" ) )

# fitted expression curves of the genes that change along pseudotime in HC (one file per lineage, made by 08_HC_Dynamics_AssociationTest_Analysis.R)
# list with: genes, raw_curve (genes x 50 bins), row_z (z-scores), tested_genes, bins
curves = list()
for( i in 1:3 )
{
    curves[[ lineageTab$unit_id[ i ] ]] = readRDS( paste0( objectPath, "HC_Dynamics_Curves_", lineageTab$trajectory_id[ i ], ".rds" ) )
}


#
### Gene modules: hierarchical clustering of the z-score curves into 4 modules per lineage
#

assignTab = c()
moduleTab = c()
for( i in 1:3 )
{
    u = lineageTab$unit_id[ i ]
    x = curves[[ u ]]$row_z

    hc = hclust( dist( x ), method = "ward.D2" )
    cluster = cutree( hc, k = 4 )

    #position of the maximum of each gene and difference between maximum and minimum of the fitted curve
    peak = max.col( x, ties.method = "first" )
    dynRange = apply( curves[[ u ]]$raw_curve, 1, function( v ) diff( range( v ) ) )

    a = data.frame( raw_cluster = as.integer( cluster ), unit_id = u, trajectory_id = lineageTab$trajectory_id[ i ], endpoint = lineageTab$endpoint[ i ], condition = "HC",
                    gene = rownames( x ), dynamic_range = as.numeric( dynRange ), gene_peak_fit_position = as.integer( peak ), gene_peak_bin = as.integer( peak ), assigned = TRUE )

    #modules are numbered by the time of their peak (C01 = earliest)
    m = c()
    for( cl in 1:4 )
    {
        genes = a$gene[ a$raw_cluster == cl ]
        m = rbind( m, data.frame( raw_cluster = cl, full_module_genes = length( genes ),
                                  mean_gene_peak = mean( a$gene_peak_fit_position[ a$raw_cluster == cl ] ),
                                  median_gene_peak = median( a$gene_peak_fit_position[ a$raw_cluster == cl ] ),
                                  median_dynamic_range = median( a$dynamic_range[ a$raw_cluster == cl ] ),
                                  module_mean_curve_peak = which.max( colMeans( x[ genes, , drop = FALSE ] ) ) ) )
    }
    m = m[ order( m$mean_gene_peak, m$median_gene_peak, m$module_mean_curve_peak, m$raw_cluster ), ]
    m$block_order = 1:4
    m$block = sprintf( "C%02d", 1:4 )
    m$unit_id = u
    m$trajectory_id = lineageTab$trajectory_id[ i ]
    m$endpoint = lineageTab$endpoint[ i ]
    m$condition = "HC"

    a$block_order = m$block_order[ match( a$raw_cluster, m$raw_cluster ) ]
    a$block = m$block[ match( a$raw_cluster, m$raw_cluster ) ]
    a = a[ order( a$block_order, a$gene_peak_fit_position, -a$dynamic_range, a$gene, method = "radix" ), ]

    #rank of the genes within their module by dynamic range (the heatmap shows at most 1000 genes per module)
    a$within_module_dynamic_rank = NA
    for( b in m$block )
    {
        a$within_module_dynamic_rank[ a$block == b ] = rank( -a$dynamic_range[ a$block == b ], ties.method = "first" )
    }
    a$full_module_row_order = 1:nrow( a )

    m$displayed_genes = sapply( m$block, function( b ) sum( a$block == b & a$within_module_dynamic_rank <= 1000 ) )
    m$display_fraction = m$displayed_genes / m$full_module_genes

    assignTab = rbind( assignTab, a )
    moduleTab = rbind( moduleTab, m[ , c( "raw_cluster", "block", "full_module_genes", "mean_gene_peak", "median_gene_peak", "median_dynamic_range", "module_mean_curve_peak",
                                          "block_order", "unit_id", "trajectory_id", "endpoint", "condition", "displayed_genes", "display_fraction" ) ] )
}

# genes shown in the heatmaps
topTab = assignTab[ assignTab$within_module_dynamic_rank <= 1000, ]
topTab$display_row_order = unlist( lapply( lineageTab$unit_id, function( u ) 1:sum( topTab$unit_id == u ) ) )

fwrite( assignTab, file = paste0( tablePath, "HC_Dynamics_Modules_Genes.txt.gz" ), sep = "\t" )
fwrite( moduleTab, file = paste0( tablePath, "HC_Dynamics_Modules_Overview.txt.gz" ), sep = "\t" )
fwrite( topTab, file = paste0( tablePath, "HC_Dynamics_Modules_GenesInHeatmap.txt.gz" ), sep = "\t" )


#
### GO biological process enrichment of each module
#

# technical genes (ribosomal and mitochondrial genes, see Data/DataFileInfo.txt; columns: gene_symbol, ensembl_id)
technicalGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/TechnicalGenes.txt" ) ) )

# GO BP gene sets of MSigDB 2026.1 (made by DataPreparation/ReferenceFiles/code/MSigDB_GeneSets.R), term id: GO id
gobp = as.data.frame( readRDS( paste0( dataPath, "ReferenceFiles/GeneSets/MSigDB_2026.1_GO_BP.rds" ) ) )
term2gene = unique( data.frame( term = gobp$gs_exact_source, gene = toupper( trimws( gobp$gene_symbol ) ) ) )
term2name = unique( data.frame( term = gobp$gs_exact_source, name = gobp$gs_name ) )
goGenes = unique( term2gene$gene )

# gene universe of a lineage: all genes that were tested for this lineage (without technical genes)
goTab = c()
goTop3 = c()
for( i in 1:3 )
{
    u = lineageTab$unit_id[ i ]
    universe = unique( toupper( trimws( curves[[ u ]]$tested_genes ) ) )
    universe = universe[ universe != "" & !universe %in% technicalGenes$gene_symbol ]

    for( b in sprintf( "C%02d", 1:4 ) )
    {
        genes = unique( toupper( trimws( assignTab$gene[ assignTab$unit_id == u & assignTab$block == b ] ) ) )
        genes = intersect( genes, universe )

        # universe: tested genes in at least one GO BP gene set
        res = enricher( gene = intersect( genes, goGenes ), universe = intersect( universe, goGenes ), TERM2GENE = term2gene, TERM2NAME = term2name,
                        pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500 )
        res = as.data.frame( res )

        tab = data.frame( unit_id = character(), trajectory_id = character(), block = character() )
        if( nrow( res ) > 0 )
        {
            tab = data.frame( unit_id = u, trajectory_id = lineageTab$trajectory_id[ i ], block = b, go_id = res$ID, term_label = tools::toTitleCase( tolower( res$Description ) ),
                              gene_ratio = res$GeneRatio, background_ratio = res$BgRatio, pvalue = res$pvalue, padj = res$p.adjust, qvalue = res$qvalue,
                              genes = gsub( "/", ";", res$geneID, fixed = TRUE ), overlap = res$Count, query_genes = length( genes ), background_genes = length( universe ) )
            tab = tab[ order( tab$padj, tab$pvalue, -tab$overlap, tab$go_id ), ]
            goTab = rbind( goTab, tab )

            #the three most significant terms are written next to the heatmap
            top3 = tab[ 1:min( 3, nrow( tab ) ), ]
            top3$display_rank = 1:nrow( top3 )
            goTop3 = rbind( goTop3, top3 )
        }
        else
        {
            goTop3 = rbind( goTop3, data.frame( unit_id = u, trajectory_id = lineageTab$trajectory_id[ i ], block = b, go_id = NA, term_label = "No enriched GO BP term",
                                                gene_ratio = NA, background_ratio = NA, pvalue = NA, padj = NA, qvalue = NA, genes = NA, overlap = NA, query_genes = NA, background_genes = NA, display_rank = NA ) )
        }
    }
}

fwrite( goTab, file = paste0( tablePath, "HC_Dynamics_Modules_GO_BP.txt.gz" ), sep = "\t" )
fwrite( goTop3, file = paste0( tablePath, "HC_Dynamics_Modules_GO_BP_Top3.txt" ), sep = "\t", na = "NA" )


#
### Distribution of the cell clusters along pseudotime (shown above the heatmaps)
#

# slingshot result (made by 01_Trajectories_Slingshot_Analysis.R)
# list with: slingshot_sce, lineages, pseudotime, cell_weights, cell_metadata
trajectory <- readRDS( paste0( objectPath, "OPC_OL_Slingshot.rds" ) )

# cells of each lineage and their bin (made by 08_HC_Dynamics_AssociationTest_Analysis.R)
# columns: trajectory_id, cell_id, pseudotime, bin_index
binCells = as.data.frame( fread( paste0( tablePath, "HC_Dynamics_Bins_Cells.txt.gz" ) ) )

#
# Smoothing of a vector with a gaussian kernel, the ends are filled with the first and last value
#
gaussianSmooth <- function( x, sigma = 1.35, radius = 4 )
{
    w = exp( -( ( -radius:radius )^2 ) / ( 2 * sigma^2 ) )
    w = w / sum( w )
    padded = c( rep( x[ 1 ], radius ), x, rep( x[ length( x ) ], radius ) )
    as.numeric( stats::filter( padded, filter = w, sides = 2 )[ 1:length( x ) + radius ] )
}

waves = list()
binTab = c()
waveTab = c()
for( i in 1:3 )
{
    u = lineageTab$unit_id[ i ]
    path = as.character( trajectory$lineages[[ i ]] )

    #HC cells of the lineage in the 50 bins of the heatmap, clusters on the path of the lineage
    d = binCells[ binCells$trajectory_id == lineageTab$trajectory_id[ i ], ]
    d$final_cluster_label = trajectory$cell_metadata$final_cluster_label[ match( d$cell_id, trajectory$cell_metadata$cell_id ) ]
    stopifnot( !anyNA( d$final_cluster_label ) )
    d = d[ d$final_cluster_label %in% path, ]
    nBins = 50

    n = table( factor( d$final_cluster_label, levels = path ), factor( d$bin_index, levels = 1:nBins ) )

    raw = c()
    smooth = c()
    for( cl in path )
    {
        nCells = as.numeric( n[ cl, ] )
        raw = rbind( raw, data.frame( bin_index = 1:nBins, final_cluster_label = cl, n_cells = nCells, cluster_peak = max( nCells ),
                                      count_scaled = if( max( nCells ) > 0 ) nCells / max( nCells ) else 0,
                                      pseudotime_scaled = ( 1:nBins - 1 ) / ( nBins - 1 ) ) )

        y = gaussianSmooth( nCells )
        yy = pmax( approx( 1:length( y ), y, xout = seq( 1, length( y ), length.out = 200 ), rule = 2 )$y, 0 )
        smooth = rbind( smooth, data.frame( final_cluster_label = cl, x_index = 1:200, pseudotime_scaled = seq( 0, 1, length.out = 200 ),
                                            n_cells_smoothed = yy, count_scaled = if( max( yy ) > 0 ) yy / max( yy ) else 0 ) )
    }
    for( tab in c( "raw", "smooth" ) )
    {
        x = get( tab )
        x$unit_id = u
        x$trajectory_id = lineageTab$trajectory_id[ i ]
        x$condition = "HC"
        x$actual_bins = nBins
        x$path_order = match( x$final_cluster_label, path )
        assign( tab, x )
    }
    binTab = rbind( binTab, raw )
    waveTab = rbind( waveTab, smooth )
    waves[[ u ]] = list( smooth = smooth, path = path )
}

fwrite( binTab, file = paste0( tablePath, "HC_Dynamics_ClusterDistribution_Bins.txt.gz" ), sep = "\t" )
fwrite( waveTab, file = paste0( tablePath, "HC_Dynamics_ClusterDistribution_Smoothed.txt.gz" ), sep = "\t" )


#
### Results for the heatmaps (10_HC_Dynamics_Heatmaps_Plots.R)
#

saveRDS( list( topTab = topTab, goTop3 = goTop3, waves = waves ), file = paste0( objectPath, "HC_Dynamics_PlotData.rds" ) )
