library(tradeSeq)
library(SingleCellExperiment)
library(Matrix)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# output folder
outputPath = paste0( home_directory, "Figure_5/" )
modulePath = paste0( outputPath, "tables/Modules/" )
fitPath = paste0( outputPath, "tables/TwoLineages/" )
tablePath = paste0( outputPath, "tables/SelectedGenes/" )
objectPath = paste0( outputPath, "objects/" )


# two genes of the respiratory chain that are in the first module of PD L1 and of PD L3: COX6B1 and ATP5F1E
disease = "PD"
genes = c( "COX6B1", "ATP5F1E" )
modelDir = "fitGAM_PD_vs_HC_TwoLineages"
# the cells shown as points are sampled (see below), every gene has its own seed
seeds = 42 + 2000 + 1:2

# results of the condition test, modules and GO enrichment (made by 02_ and 03_ of this folder)
condTab <- readRDS( paste0( objectPath, "Disease_Dynamics_ConditionTest.rds" ) )
moduleGenes = as.data.frame( fread( paste0( modulePath, "Disease_Dynamics_Modules_Genes.txt.gz" ) ) )
goTab = as.data.frame( fread( paste0( modulePath, "Disease_Dynamics_Modules_GO_BP.txt.gz" ) ) )

# test results of the two genes and are they part of the OXPHOS term of the first modules?
statTab = condTab[ condTab$gene %in% genes & condTab$disease == disease, c( "gene", "unit_id", "disease", "wald_stat", "df", "pvalue", "padj_BH", "XDE_FDR05" ) ]
statTab = statTab[ order( statTab$gene, statTab$unit_id ), ]
statTab$module_id = moduleGenes$module_id[ match( paste( statTab$unit_id, statTab$gene ), paste( moduleGenes$unit_id, moduleGenes$gene ) ) ]
oxphosTerm = goTab[ goTab$term_name == "GOBP_OXIDATIVE_PHOSPHORYLATION", ]
statTab$GO_BP_OXPHOS_FDR = oxphosTerm$padj_within_module_BH[ match( statTab$module_id, oxphosTerm$module_id ) ]
statTab$GO_BP_overlap_gene = sapply( 1:nrow( statTab ), function( k ) statTab$gene[ k ] %in% strsplit( oxphosTerm$overlap_genes[ match( statTab$module_id[ k ], oxphosTerm$module_id ) ], ";", fixed = TRUE )[[ 1 ]] )

fwrite( statTab, file = paste0( tablePath, "SelectedGenes_OXPHOS_", disease, "_Statistics.txt" ), sep = "\t" )

# inputs of the models (made by 01_TwoLineages_fitGAM_Analysis.R)
# list with: counts, pseudotime, cell_weights, cell_metadata, conditions, offset
inputs <- readRDS( paste0( objectPath, "fitGAM_Inputs_", disease, "_vs_HC.rds" ) )
# columns: gene, gene_rank, detected_cells, total_counts, batch_id
batchTab = as.data.frame( fread( paste0( fitPath, "fitGAM_GeneBatches_", disease, "_vs_HC.txt.gz" ) ) )
# compared pseudotime range of each lineage
rangeTab = as.data.frame( fread( paste0( fitPath, "Disease_Dynamics_ComparedPseudotimeRange.txt" ) ) )
rangeTab <- rangeTab[ rangeTab$disease == disease, ]


#
### Fitted curves and expression of the single cells within the compared pseudotime range
#

curveTab = c()
pointTab = c()
samplingTab = c()
for( k in 1:length( genes ) )
{
    g = genes[ k ]
    b = batchTab$batch_id[ batchTab$gene == g ]
    fit = readRDS( paste0( objectPath, modelDir, "/fit_sce_batch_", sprintf( "%04d", b ), ".rds" ) )

    #fitted curve: 1000 points per lineage and condition
    smooth = as.data.frame( predictSmooth( fit, gene = g, nPoints = 1000, tidy = TRUE ) )
    smooth$lineage_id = c( "L1", "L3" )[ as.integer( smooth$lineage ) ]
    for( i in 1:2 )
    {
        s = smooth[ smooth$lineage_id == rangeTab$lineage_id[ i ] & smooth$time >= rangeTab$test_time_start[ i ] & smooth$time <= rangeTab$test_time_end[ i ], ]
        curveTab = rbind( curveTab, data.frame( pair = paste0( disease, "_vs_HC" ), disease = disease, gene = g, lineage_id = s$lineage_id, condition = as.character( s$condition ),
                                                condition_plot = ifelse( s$condition == "HC", "HC", "Disease" ), pseudotime = s$time, fitted_log_expression = log1p( s$yhat ) ) )
    }

    #cells: design matrix of the model tells to which lineage a cell was assigned (columns l1_1, l1_2, l2_1, l2_2) and its pseudotime (t1, t2)
    dm = as.data.frame( colData( fit )$tradeSeq$dm )
    geneCounts = as.numeric( inputs$counts[ g, ] )
    cellPoints = c()
    for( i in 1:2 )
    {
        for( ci in 1:2 )
        {
            cond = levels( inputs$conditions )[ ci ]
            keep = dm[ , paste0( "l", i, "_", ci ) ] == 1 & inputs$cell_metadata$final_cluster_label != "OPC" &
                   dm[ , paste0( "t", i ) ] >= rangeTab$test_time_start[ i ] & dm[ , paste0( "t", i ) ] <= rangeTab$test_time_end[ i ]
            cellPoints = rbind( cellPoints, data.frame( pair = paste0( disease, "_vs_HC" ), disease = disease, gene = g, cell_id = colnames( inputs$counts )[ keep ], lineage_id = rangeTab$lineage_id[ i ],
                                                        condition = cond, condition_plot = ifelse( cond == "HC", "HC", "Disease" ),
                                                        pseudotime = dm[ keep, paste0( "t", i ) ], log_expression = log1p( geneCounts[ keep ] ) ) )
        }
    }

    #half of the cells are shown (randomly selected), otherwise the curves are hidden behind the points
    set.seed( seeds[ k ] )
    shown = cellPoints[ sample.int( nrow( cellPoints ), size = floor( nrow( cellPoints ) * 0.5 ), replace = FALSE ), ]
    pointTab = rbind( pointTab, shown )
    samplingTab = rbind( samplingTab, data.frame( gene = g, eligible_point_rows = nrow( cellPoints ), sampled_point_rows = nrow( shown ), sampling_fraction_realized = nrow( shown ) / nrow( cellPoints ),
                                                  unique_eligible_cells = length( unique( cellPoints$cell_id ) ), unique_sampled_cells = length( unique( shown$cell_id ) ), fit_batch_id = b ) )
}

fwrite( samplingTab, file = paste0( tablePath, "SelectedGenes_OXPHOS_", disease, "_ShownCells.txt" ), sep = "\t" )

# fitted curves and cells for the plot (06_SelectedGenes_OXPHOS_PD_Plots.R)
saveRDS( list( curveTab = curveTab, pointTab = pointTab ), file = paste0( objectPath, "SelectedGenes_OXPHOS_", disease, "_PlotData.rds" ) )
