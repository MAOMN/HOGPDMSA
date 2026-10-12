library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# transcription factor ranking of our PD comparison (Figure 3G)
patientPath = paste0( home_directory, "Figure_3/" )

# output folder
outputPath = paste0( home_directory, "SupplementaryTables/Table_S1/" )
tablePath = paste0( outputPath, "tables/" )


#
### Selected TFs: for each data set its own 10 most consistently significant TFs plus the 10 TFs shown for our PD comparison (Figure 3G)
#

cohorts = c( martirosyan = "Martirosyan", smajic = "Smajic", lee = "Lee" )

# Supplementary Table S1 (made by 01_TF_Activity_ExternalPDCohorts_Analysis.R), read as text, the fields of the TFs that were not tested are empty
tableS1 = as.data.frame( fread( paste0( tablePath, "Supplementary_Table_S1_TF_activity.txt" ), colClasses = "character", na.strings = NULL ) )
tableS1$FDR = suppressWarnings( as.numeric( tableS1$FDR_BH_within_population ) )

# the 10 TFs of Figure 3G (made by 14_TF_Activity_OPC_OL_Analysis.R of Figure 3)
# columns: source, FDR_sum, eligible, missing_groups, FDR_OL1, ..., display_rank
ownRanking = as.data.frame( fread( paste0( patientPath, "tables/TF_Activity/TF_Activity_OPC_OL_Ranking_PD_FDRSum.txt" ) ) )
ownTop = ownRanking$source[ 1:10 ]

selectionTab = c()
for( ds in cohorts )
{
    x = tableS1[ tableS1$dataset == ds, ]

    #ranking: sum of the FDR values of all populations of the data set, TFs that could not be tested in all populations are not ranked
    tfs = unique( x$TF )
    fdrSum = sapply( tfs, function( tf ) sum( x$FDR[ x$TF == tf ] ) )
    ranking = data.frame( TF = tfs, FDR_sum = fdrSum )
    ranking = ranking[ is.finite( ranking$FDR_sum ), ]
    ranking = ranking[ order( ranking$FDR_sum, ranking$TF, method = "radix" ), ]
    cohortTop = ranking$TF[ 1:10 ]
    selected = c( cohortTop, ownTop[ !ownTop %in% cohortTop ] )

    selectionTab = rbind( selectionTab, data.frame( dataset = ds, TF = selected, display_order = 1:length( selected ), in_dataset_top10 = selected %in% cohortTop, in_NatComm_PD_top10 = selected %in% ownTop,
                                                    dataset_FDR_sum = ranking$FDR_sum[ match( selected, ranking$TF ) ] ) )
}

fwrite( selectionTab, file = paste0( tablePath, "ExternalPDCohorts_SelectedTFs.txt" ), sep = "\t", quote = FALSE )
