library(speckle)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_3/" )
tablePath = paste0( outputPath, "tables/" )
objectPath = paste0( outputPath, "objects/" )


# table of the clustered OPCs and oligodendrocytes (03_OPC_OL_Clustering_Analysis.R), only cells of the kept clusters
cells = as.data.frame( fread( paste0( tablePath, "OPC_OL_Clustering/OPC_OL_CellTable.txt.gz" ) ) )
cells <- cells[ cells$final_retained, ]

# donor information (columns: sample_id, Diagnosis_author, Diagnosis, Sex, Age, PMI, brain_bank, extract, flowcell, disease_duration, subtype)
# made by DataPreparation/snRNAseq_HC_MSA_PD/code/DonorCovariates_AuthorMetadata.R from the data of Rydbirk et al. (2026)
donor_info = as.data.frame( fread( paste0( dataPath, "snRNAseq_HC_MSA_PD/NatComm_DonorCovariates.txt" ) ) )

states = c( "OPC", paste0( "OL", 1:6 ) )
diagnoses = c( "HC", "MSA", "PD" )


#
### Proportion of each cell state per donor
#

donors = unique( cells[ , c( "diagnosis", "sample_id" ) ] )
donors$diagnosis <- factor( donors$diagnosis, levels = diagnoses )
donors <- donors[ order( donors$diagnosis, donors$sample_id, method = "radix" ), ]

cellNumbers = table( cells$sample_id, cells$final_cluster_label )

donorProp = c()
for( i in 1:nrow( donors ) )
{
    n = as.numeric( cellNumbers[ donors$sample_id[ i ], states ] )
    donorProp = rbind( donorProp, data.frame( diagnosis = donors$diagnosis[ i ], sample_id = donors$sample_id[ i ], final_cluster_label = states,
                                              total_retained_cells = sum( n ), n_cells = n, donor_proportion = n / sum( n ) ) )
}
donorProp$final_cluster_label <- factor( donorProp$final_cluster_label, levels = states )

fwrite( donorProp, file = paste0( tablePath, "CellStateProportions/OPC_OL_DonorProportions.txt" ), sep = "\t" )

# mean proportion over the donors of each diagnosis
composition = c()
for( dx in diagnoses )
{
    for( s in states )
    {
        d = donorProp[ donorProp$diagnosis == dx & donorProp$final_cluster_label == s, ]
        composition = rbind( composition, data.frame( diagnosis = dx, final_cluster_label = s, mean_prop = mean( d$donor_proportion ),
                                                      n_donors = nrow( d ), pooled_cells = sum( d$n_cells ) ) )
    }
}
composition$diagnosis <- factor( composition$diagnosis, levels = diagnoses )
composition$final_cluster_label <- factor( composition$final_cluster_label, levels = states )
composition$group_total_cells <- sapply( composition$diagnosis, function( dx ) sum( composition$pooled_cells[ composition$diagnosis == dx ] ) )
composition$pooled_prop <- composition$pooled_cells / composition$group_total_cells


#
### Test of the proportions with propeller: each disease against HC, adjusted for age
#

statistics = c()
for( disease in c( "MSA", "PD" ) )
{
    d = cells[ cells$diagnosis %in% c( "HC", disease ), ]
    ids = sort( unique( d$sample_id ) )

    info = donor_info[ match( ids, donor_info$sample_id ), ]
    info$Diagnosis <- factor( info$Diagnosis, levels = c( "HC", disease ) )
    info$Age_z <- as.numeric( scale( info$Age ) )

    design <- model.matrix( ~ 0 + Diagnosis + Age_z, data = info )
    rownames( design ) <- ids
    # disease minus HC
    contrast <- c( -1, 1, 0 )
    names( contrast ) <- colnames( design )

    # logit transformed proportions of each donor
    props <- getTransformedProps( clusters = d$final_cluster_label, sample = d$sample_id, transform = "logit" )

    res <- propeller.ttest( prop.list = props, design = design, contrasts = contrast, robust = TRUE, trend = FALSE, sort = FALSE )

    # save the proportions of the donors and the test result as R object
    saveRDS( list( props = props, donor_info = info, result = res ), file = paste0( objectPath, "Propeller_", disease, "_vs_HC.rds" ) )
    fwrite( data.frame( state = rownames( res ), res ), file = paste0( tablePath, "CellStateProportions/Propeller_", disease, "_vs_HC_Result.txt" ), sep = "\t" )

    res <- res[ states, ]
    statistics = rbind( statistics, data.frame( state = states,
                                                mean_HC = res$PropMean.DiagnosisHC,
                                                mean_disease = res[ , paste0( "PropMean.Diagnosis", disease ) ],
                                                PropRatio = res$PropRatio, Tstatistic = res$Tstatistic, P.Value = res$P.Value, FDR = res$FDR,
                                                comparison = paste0( disease, "_vs_HC" ), diagnosis = disease,
                                                n_HC = sum( info$Diagnosis == "HC" ), n_disease = sum( info$Diagnosis == disease ) ) )
}
statistics$raw_difference_percentage_points <- 100 * ( statistics$mean_disease - statistics$mean_HC )
statistics$significance <- sapply( statistics$FDR, function( x ) if ( x < 0.001 ) "***" else if ( x < 0.01 ) "**" else if ( x < 0.05 ) "*" else "" )

fwrite( statistics, file = paste0( tablePath, "CellStateProportions/Propeller_Statistics.txt" ), sep = "\t" )

# data of the plot in 06_CellStateProportions_Plots.R
saveRDS( list( composition = composition, statistics = statistics ), file = paste0( objectPath, "Propeller_PlotData.rds" ) )
