library(tradeSeq)
library(SingleCellExperiment)
library(data.table)
library(parallel)

# change to your project directory
home_directory <- "path/to/project_directory/"

# output folder
outputPath = paste0( home_directory, "Figure_4/" )
fitPath = paste0( outputPath, "tables/fitGAM/" )
tablePath = paste0( outputPath, "tables/HC_Dynamics/" )
objectPath = paste0( outputPath, "objects/" )


# which genes change along pseudotime in the healthy controls?
# the fitted models (three lineages, conditions HC, MSA and PD) are stored in batches of 100 genes (made by 07_fitGAM_Analysis.R)

# columns: gene, batch_id (made by 07_fitGAM_Analysis.R)
batchTab = as.data.frame( fread( paste0( fitPath, "fitGAM_GeneBatches_ThreeLineages.txt.gz" ) ) )
batches = unique( batchTab$batch_id )

fitFile <- function( b )
{
    paste0( objectPath, "fitGAM_HC_MSA_PD_ThreeLineages/fit_sce_batch_", sprintf( "%04d", b ), ".rds" )
}

lineageTab = data.frame( unit_id = c( "HC_L1_TDE", "HC_L2_TDE", "HC_L3_TDE" ), trajectory_id = c( "L1", "L2", "L3" ), model_lineage = 1:3, endpoint = c( "OL5", "OL6", "OL3" ) )


#
### Association test of tradeSeq
#

# tested against a log2 fold change of 1 along the lineage, separately for each lineage and condition
testBatch <- function( b )
{
    fit = readRDS( fitFile( b ) )
    genes = as.character( rowData( fit )$tradeSeq$name )
    stopifnot( identical( genes, batchTab$gene[ batchTab$batch_id == b ] ) )
    res = associationTest( fit, lineages = TRUE, l2fc = 1, inverse = "eigen" )
    data.frame( batch_id = b, gene = genes, converged = rowData( fit )$tradeSeq$converged, res, row.names = NULL )
}
assocTab = do.call( rbind, mclapply( batches, testBatch, mc.cores = 8, mc.preschedule = FALSE ) )

fwrite( assocTab, file = paste0( tablePath, "HC_Dynamics_AssociationTest_AllLineagesAndConditions.txt.gz" ), sep = "\t" )

# results of the HC condition for the three lineages, p-values adjusted within each lineage
hcTab = c()
overview = c()
for( i in 1:3 )
{
    tab = data.frame( assocTab[ , c( "batch_id", "gene", "converged" ) ],
                      waldStat = assocTab[ , paste0( "waldStat_lineage", i, "_conditionHC" ) ],
                      df = assocTab[ , paste0( "df_lineage", i, "_conditionHC" ) ],
                      pvalue = assocTab[ , paste0( "pvalue_lineage", i, "_conditionHC" ) ] )
    # genes where the model did not converge or without p-value are not used
    tab$eligible_for_BH = tab$converged %in% TRUE & is.finite( tab$pvalue )
    tab$padj_BH_lineage = NA
    tab$padj_BH_lineage[ tab$eligible_for_BH ] = p.adjust( tab$pvalue[ tab$eligible_for_BH ], method = "BH" )
    tab$significant_TDE = tab$eligible_for_BH & is.finite( tab$padj_BH_lineage ) & tab$padj_BH_lineage < 0.05
    tab$unit_id = lineageTab$unit_id[ i ]
    tab$trajectory_id = lineageTab$trajectory_id[ i ]
    tab$model_lineage = i
    tab$endpoint = lineageTab$endpoint[ i ]
    tab$condition = "HC"
    hcTab = rbind( hcTab, tab )

    overview = rbind( overview, data.frame( unit_id = lineageTab$unit_id[ i ], trajectory_id = lineageTab$trajectory_id[ i ], endpoint = lineageTab$endpoint[ i ],
                                            candidate_genes = nrow( tab ), converged_genes = sum( tab$converged %in% TRUE ), BH_family_genes = sum( tab$eligible_for_BH ),
                                            significant_TDE_genes = sum( tab$significant_TDE ),
                                            minimum_pvalue = min( tab$pvalue[ tab$eligible_for_BH ] ), minimum_padj = min( tab$padj_BH_lineage[ tab$eligible_for_BH ] ), FDR_cutoff = 0.05 ) )
}

fwrite( hcTab, file = paste0( tablePath, "HC_Dynamics_AssociationTest_HC.txt.gz" ), sep = "\t" )
fwrite( overview, file = paste0( tablePath, "HC_Dynamics_NumberOfDynamicGenes.txt" ), sep = "\t" )


#
### Fitted expression curves of the significant genes (HC) in 50 bins with the same number of cells
#
# - cells of a lineage: the HC cells that the model assigns to the lineage (design matrix column l<lineage>_1 = 1), ranked by pseudotime
# - the curve of each bin is taken at the median pseudotime of its cells (predictSmooth with 1000 points, linear interpolation)
# - z-score of each gene over the 50 bins

sigGenes = lapply( lineageTab$unit_id, function( u ) hcTab$gene[ hcTab$unit_id == u & hcTab$significant_TDE ] )
nBins = 50

# cells of the trajectories, in the order of the cells of the models (made by 01_Trajectories_Slingshot_Analysis.R)
trajectory <- readRDS( paste0( objectPath, "OPC_OL_Slingshot.rds" ) )
fit1 = readRDS( fitFile( 1 ) )
dm = colData( fit1 )$tradeSeq$dm
stopifnot( levels( fit1$tradeSeq$conditions )[ 1 ] == "HC", nrow( dm ) == nrow( trajectory$cell_metadata ) )
cellTab = c()
binTab = c()
binTime = list()
for( i in 1:3 )
{
    k = dm[ , paste0( "l", i, "_1" ) ] == 1
    t = dm[ k, paste0( "t", i ) ]
    bin = ceiling( rank( t, ties.method = "first" ) / length( t ) * nBins )
    binTime[[ i ]] = as.numeric( tapply( t, bin, median ) )
    cellTab = rbind( cellTab, data.frame( trajectory_id = lineageTab$trajectory_id[ i ], cell_id = trajectory$cell_metadata$cell_id[ k ], pseudotime = t, bin_index = bin ) )
    binTab = rbind( binTab, data.frame( trajectory_id = lineageTab$trajectory_id[ i ], bin_index = 1:nBins, HC_cells = as.integer( table( factor( bin, levels = 1:nBins ) ) ),
                                        pseudotime_min = as.numeric( tapply( t, bin, min ) ), pseudotime_max = as.numeric( tapply( t, bin, max ) ), pseudotime_median = binTime[[ i ]] ) )
}
fwrite( cellTab, file = paste0( tablePath, "HC_Dynamics_Bins_Cells.txt.gz" ), sep = "\t" )
fwrite( binTab, file = paste0( tablePath, "HC_Dynamics_Bins.txt" ), sep = "\t" )

curvesOfBatch <- function( b )
{
    wanted = lapply( sigGenes, intersect, batchTab$gene[ batchTab$batch_id == b ] )
    genes = unique( unlist( wanted ) )
    res = list()
    if( length( genes ) > 0 )
    {
        fit = readRDS( fitFile( b ) )
        smooth = as.data.frame( predictSmooth( fit, gene = genes, nPoints = 1000, tidy = TRUE ) )
        smooth$gene = as.character( smooth$gene )
        for( i in 1:3 )
        {
            g = wanted[[ i ]]
            if( length( g ) > 0 )
            {
                d = smooth[ smooth$condition == "HC" & smooth$lineage == i & smooth$gene %in% g, ]
                mat = t( sapply( g, function( x ) {
                    e = d[ d$gene == x, ]
                    e = e[ order( e$time ), ]
                    stopifnot( nrow( e ) == 1000, min( e$time ) <= min( binTime[[ i ]] ), max( e$time ) >= max( binTime[[ i ]] ) )
                    approx( e$time, e$yhat, xout = binTime[[ i ]] )$y } ) )
                colnames( mat ) = sprintf( "bin_%02d", 1:nBins )
                res[[ lineageTab$unit_id[ i ] ]] = mat
            }
        }
    }
    res
}
curveList = mclapply( batches, curvesOfBatch, mc.cores = 8, mc.preschedule = FALSE )
stopifnot( !any( sapply( curveList, inherits, "try-error" ) ) )

for( i in 1:3 )
{
    u = lineageTab$unit_id[ i ]
    raw = do.call( rbind, lapply( curveList, function( x ) x[[ u ]] ) )
    raw = raw[ sigGenes[[ i ]], ]

    # z-score of each gene over the 50 bins (genes with constant curve are removed)
    geneMean = rowMeans( raw )
    geneSd = apply( raw, 1, sd )
    valid = is.finite( geneMean ) & is.finite( geneSd ) & geneSd > sqrt( .Machine$double.eps )
    zscore = ( raw - geneMean ) / geneSd

    saveRDS( list( unit_id = u, trajectory_id = lineageTab$trajectory_id[ i ], model_lineage = i, endpoint = lineageTab$endpoint[ i ],
                   genes = rownames( raw )[ valid ], raw_curve = raw[ valid, ], row_z = zscore[ valid, ],
                   tested_genes = hcTab$gene[ hcTab$unit_id == u & hcTab$eligible_for_BH ],
                   bins = binTab[ binTab$trajectory_id == lineageTab$trajectory_id[ i ], ] ),
             file = paste0( objectPath, "HC_Dynamics_Curves_", lineageTab$trajectory_id[ i ], ".rds" ) )
}
