library(tradeSeq)
library(SingleCellExperiment)
library(data.table)
library(parallel)

# change to your project directory
home_directory <- "path/to/project_directory/"

# output folder
outputPath = paste0( home_directory, "Figure_5/" )
fitPath = paste0( outputPath, "tables/TwoLineages/" )
tablePath = paste0( outputPath, "tables/ConditionTest/" )
objectPath = paste0( outputPath, "objects/" )


# which genes have a different expression course in disease compared to HC?
# the fitted models (two lineages, conditions HC and disease) are stored in batches of 100 genes (made by 01_TwoLineages_fitGAM_Analysis.R)

modelDir = c( MSA = "fitGAM_MSA_vs_HC_TwoLineages", PD = "fitGAM_PD_vs_HC_TwoLineages" )
fitFile <- function( disease, b )
{
    paste0( objectPath, modelDir[ disease ], "/fit_sce_batch_", sprintf( "%04d", b ), ".rds" )
}

# pseudotime range (knots) that is compared for each lineage (made by 01_TwoLineages_fitGAM_Analysis.R)
# columns: disease, model_lineage_index, original_lineage_index, lineage_id, endpoint, ..., test_knot_start, test_knot_end, test_time_start, test_time_end
rangeTab <- readRDS( paste0( objectPath, "Disease_Dynamics_ComparedPseudotimeRange.rds" ) )
rangeTab$unit_id <- paste0( rangeTab$disease, "_", rangeTab$lineage_id )

batchTabs = list()
for( disease in c( "MSA", "PD" ) )
{
    # columns: gene, gene_rank, detected_cells, total_counts, batch_id (made by 01_TwoLineages_fitGAM_Analysis.R)
    batchTabs[[ disease ]] = as.data.frame( fread( paste0( fitPath, "fitGAM_GeneBatches_", disease, "_vs_HC.txt.gz" ) ) )
}


#
### Condition test of tradeSeq
#

# tested against a log2 fold change of 1, separately for each lineage, between the knots given above
tasks = rbind( data.frame( disease = "MSA", batch_id = unique( batchTabs$MSA$batch_id ) ), data.frame( disease = "PD", batch_id = unique( batchTabs$PD$batch_id ) ) )

testBatch <- function( k )
{
    disease = tasks$disease[ k ]
    b = tasks$batch_id[ k ]
    fit = readRDS( fitFile( disease, b ) )
    genes = batchTabs[[ disease ]]$gene[ batchTabs[[ disease ]]$batch_id == b ]
    stopifnot( identical( rownames( fit ), genes ) )

    res = c()
    for( i in 1:2 )
    {
        r = rangeTab[ rangeTab$disease == disease & rangeTab$model_lineage_index == i, ]
        test = conditionTest( fit, global = FALSE, lineages = TRUE, knots = c( r$test_knot_start, r$test_knot_end ), l2fc = 1 )
        res = rbind( res, data.frame( gene = genes, batch_id = b, converged = as.logical( rowData( fit )$tradeSeq$converged ),
                                      wald_stat = test[ , paste0( "waldStat_lineage", i ) ], df = test[ , paste0( "df_lineage", i ) ], pvalue = test[ , paste0( "pvalue_lineage", i ) ],
                                      unit_id = r$unit_id ) )
    }
    res
}
testTab = do.call( rbind, mclapply( 1:nrow( tasks ), testBatch, mc.cores = 4, mc.preschedule = FALSE ) )

# p-values adjusted within each disease and lineage
condTab = c()
for( k in 1:nrow( rangeTab ) )
{
    r = rangeTab[ k, ]
    tab = testTab[ testTab$unit_id == r$unit_id, c( "gene", "batch_id", "converged", "wald_stat", "df", "pvalue" ) ]
    tab = tab[ order( tab$batch_id ), ]
    # genes where the model did not converge or without valid test are not used
    tab$finite_test = tab$converged %in% TRUE & is.finite( tab$wald_stat ) & is.finite( tab$df ) & is.finite( tab$pvalue ) & tab$pvalue >= 0 & tab$pvalue <= 1
    tab$padj_BH = NA
    tab$padj_BH[ tab$finite_test ] = p.adjust( tab$pvalue[ tab$finite_test ], method = "BH" )
    tab$unit_id = r$unit_id
    tab$disease = r$disease
    tab$reference = "HC"
    tab$lineage_id = r$lineage_id
    tab$endpoint = r$endpoint
    tab$model_lineage_index = r$model_lineage_index
    tab$original_lineage_index = r$original_lineage_index
    tab$test_knot_start = r$test_knot_start
    tab$test_knot_end = r$test_knot_end
    tab$XDE_FDR05 = is.finite( tab$padj_BH ) & tab$padj_BH < 0.05
    condTab = rbind( condTab, tab )
}

overview = data.frame( unit_id = rangeTab$unit_id,
                       genes = sapply( rangeTab$unit_id, function( u ) sum( condTab$unit_id == u ) ),
                       finite_tests = sapply( rangeTab$unit_id, function( u ) sum( condTab$finite_test[ condTab$unit_id == u ] ) ),
                       XDE_FDR05 = sapply( rangeTab$unit_id, function( u ) sum( condTab$XDE_FDR05[ condTab$unit_id == u ] ) ) )

fwrite( condTab, file = paste0( tablePath, "Disease_Dynamics_ConditionTest.txt.gz" ), sep = "\t" )
fwrite( overview, file = paste0( tablePath, "Disease_Dynamics_NumberOfDifferentGenes.txt" ), sep = "\t" )
saveRDS( condTab, file = paste0( objectPath, "Disease_Dynamics_ConditionTest.rds" ) )


#
### Fitted expression curves of the significant genes: HC and disease within the compared pseudotime range
#
# - cells of a lineage: the cells that the model assigns to the lineage (HC and disease, design matrix columns l<i>_1 and l<i>_2) within the compared range, ranked by pseudotime
# - 50 groups with the same number of cells; the curves of HC and disease are taken at the median pseudotime of each group (predictSmooth with 1000 points, linear interpolation)
# - log1p, z-score of each gene over its 100 values (50 HC, 50 disease)

nGroups = 50
groupTab = c()
for( disease in c( "MSA", "PD" ) )
{
    # the design matrix is the same in all batches of a model
    dm = colData( readRDS( fitFile( disease, 1 ) ) )$tradeSeq$dm
    for( k in which( rangeTab$disease == disease ) )
    {
        r = rangeTab[ k, ]
        i = r$model_lineage_index
        inLineage = dm[ , paste0( "l", i, "_1" ) ] == 1 | dm[ , paste0( "l", i, "_2" ) ] == 1
        t = dm[ inLineage, paste0( "t", i ) ]
        t = sort( t[ t >= r$test_time_start & t <= r$test_time_end ] )
        group = ceiling( seq_along( t ) / length( t ) * nGroups )
        groupTab = rbind( groupTab, data.frame( unit_id = r$unit_id, group = 1:nGroups, cells = as.integer( table( group ) ),
                                                pseudotime_from = as.numeric( tapply( t, group, min ) ), pseudotime_to = as.numeric( tapply( t, group, max ) ),
                                                median_pseudotime = as.numeric( tapply( t, group, median ) ) ) )
    }
}
groupTab$position = ( groupTab$median_pseudotime - rangeTab$test_time_start[ match( groupTab$unit_id, rangeTab$unit_id ) ] ) /
                    ( rangeTab$test_time_end - rangeTab$test_time_start )[ match( groupTab$unit_id, rangeTab$unit_id ) ]
fwrite( groupTab, file = paste0( tablePath, "Disease_Dynamics_CellGroups.txt" ), sep = "\t" )

sigTab = condTab[ condTab$XDE_FDR05, ]
curveTasks = unique( sigTab[ , c( "unit_id", "batch_id" ) ] )

curvesOfBatch <- function( k )
{
    u = curveTasks$unit_id[ k ]
    b = curveTasks$batch_id[ k ]
    r = rangeTab[ rangeTab$unit_id == u, ]
    genes = sigTab$gene[ sigTab$unit_id == u & sigTab$batch_id == b ]

    fit = readRDS( fitFile( r$disease, b ) )
    smooth = as.data.frame( predictSmooth( fit, gene = genes, nPoints = 1000, tidy = TRUE ) )
    smooth = smooth[ smooth$lineage == r$model_lineage_index, ]

    targetTime = groupTab$median_pseudotime[ groupTab$unit_id == u ]
    res = list()
    for( cond in c( "HC", r$disease ) )
    {
        mat = t( sapply( genes, function( g ) {
          d = smooth[ smooth$gene == g & smooth$condition == cond, ]
          d = d[ order( d$time ), ]
          stopifnot( min( d$time ) <= min( targetTime ), max( d$time ) >= max( targetTime ) )
          approx( x = d$time, y = d$yhat, xout = targetTime, method = "linear", ties = "ordered" )$y } ) )
        rownames( mat ) = genes
        res[[ if( cond == "HC" ) "HC" else "disease" ]] = mat
    }
    res
}
curveList = mclapply( 1:nrow( curveTasks ), curvesOfBatch, mc.cores = 4, mc.preschedule = FALSE )
stopifnot( !any( sapply( curveList, inherits, "try-error" ) ) )

curves = list()
for( k in 1:nrow( rangeTab ) )
{
    r = rangeTab[ k, ]
    u = r$unit_id
    genes = sigTab$gene[ sigTab$unit_id == u ]
    parts = curveList[ which( curveTasks$unit_id == u ) ]
    parts = parts[ order( curveTasks$batch_id[ curveTasks$unit_id == u ] ) ]

    hc = do.call( rbind, lapply( parts, function( x ) x$HC ) )[ genes, ]
    dis = do.call( rbind, lapply( parts, function( x ) x$disease ) )[ genes, ]
    raw = cbind( hc, dis )
    colnames( raw ) = c( sprintf( "HC_%03d", 1:nGroups ), sprintf( "%s_%03d", r$disease, 1:nGroups ) )

    # log expression, z-score of each gene (HC and disease curve together)
    logExpr = log1p( raw )
    zscore = ( logExpr - rowMeans( logExpr ) ) / apply( logExpr, 1, sd )
    stopifnot( all( is.finite( zscore ) ) )

    g = groupTab[ groupTab$unit_id == u, ]
    curves[[ u ]] = list( unit_id = u, disease = r$disease, lineage_id = r$lineage_id, endpoint = r$endpoint, genes = genes,
                          raw_response = raw, log1p_response = logExpr, scaled = zscore,
                          time = g$position, median_pseudotime = g$median_pseudotime )
}

saveRDS( curves, file = paste0( objectPath, "Disease_Dynamics_Curves.rds" ) )
