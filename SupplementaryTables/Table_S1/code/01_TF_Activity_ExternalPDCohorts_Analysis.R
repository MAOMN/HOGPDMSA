library(DESeq2)
library(decoupleR)
library(viper)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/" )

# network (processed for the HOG cells, made by 05_TF_Activity_Analysis.R of Figure 2)
hogPath = paste0( home_directory, "Figure_2/" )

# cell populations of the published data sets (made by the scripts of DataPreparation/ExternalPDCohorts)
cohortPath = paste0( home_directory, "DataPreparation/ExternalPDCohorts/tables/" )

# output folder
outputPath = paste0( home_directory, "SupplementaryTables/Table_S1/" )
tablePath = paste0( outputPath, "tables/" )


# Supplementary Table S1: transcription factor activities (PD vs HC) in three published snRNA-seq data sets
# same analysis as for our own cell states (Figure 3G): pseudobulk per donor, DESeq2 with age in the model, VIPER with the brainSCOPE oligodendrocyte network

cohorts = c( martirosyan = "Martirosyan", smajic = "Smajic", lee = "Lee" )

# cell populations of each data set: OPC and the oligodendrocyte clusters (OL1 = largest)
# columns: population, sample_id, diagnosis, sex, age, n_cells
populations = list()
for( cohort in names( cohorts ) )
{
    donorTab = as.data.frame( fread( paste0( cohortPath, cohorts[ cohort ], "_OPC_OL_CellsPerDonorAndPopulation.txt" ) ) )
    populations[[ cohort ]] = c( "OPC", sort( setdiff( unique( donorTab$population ), "OPC" ) ) )
}

# technical genes (ribosomal and mitochondrial genes, see Data/DataFileInfo.txt; columns: gene_symbol, ensembl_id, genes without symbol are named by their Ensembl id in the Smajic data set)
technicalGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/TechnicalGenes.txt" ) ) )
technicalNames = c( technicalGenes$gene_symbol, technicalGenes$ensembl_id[ technicalGenes$ensembl_id != "" ] )

# signed brainSCOPE oligodendrocyte network
# columns: source, target, mor, likelihood
net = as.data.frame( fread( paste0( hogPath, "tables/TF_Activity/brainSCOPE_Oli_Network_SignedLinks.txt.gz" ) ) )
targetsPerTF = split( net$target, net$source )


#
### Differential expression PD vs HC and transcription factor activities for each data set and cell population
#

tfTable = c()
designTab = c()
for( cohort in names( cohorts ) )
{
    for( population in populations[[ cohort ]] )
    {
        id = paste( cohort, "PD_vs_HC", population, sep = "__" )

        # pseudobulk counts per donor of one cell population (cell populations as defined within each data set)
        # list with: counts (genes x donors), coldata (sample_id, diagnosis, sex, age, n_cells)
        pb <- readRDS( paste0( dataPath, "ExternalPDCohorts/Pseudobulk/", cohort, "/", population, ".rds" ) )

        # donors with at least 10 cells
        d = as.data.frame( pb$coldata )
        d = d[ d$n_cells >= 10 & d$diagnosis %in% c( "HC", "PD" ), ]
        d = d[ order( d$sample_id, method = "radix" ), ]
        d$Diagnosis = factor( d$diagnosis, levels = c( "HC", "PD" ) )
        d$Age_z = as.numeric( scale( as.numeric( d$age ) ) )
        rownames( d ) = d$sample_id

        counts = pb$counts[ , d$sample_id ]
        # gene has to be detected in at least half of the donors
        keep = rowSums( counts > 0 ) >= ceiling( ncol( counts ) / 2 ) & rowSums( counts ) >= 10

        dds <- DESeqDataSetFromMatrix( countData = round( counts[ keep, ] ), colData = d[ , c( "sample_id", "Age_z", "Diagnosis" ) ], design = ~ Age_z + Diagnosis )
        dds <- DESeq( dds, fitType = "local", quiet = TRUE )
        res <- as.data.frame( results( dds, contrast = c( "Diagnosis", "PD", "HC" ), alpha = 0.05 ) )

        resTable = data.frame( gene = rownames( res ), res, row.names = NULL )
        # technical genes are not used
        resTable$technical_excluded = toupper( resTable$gene ) %in% technicalNames
        resTable$TF_input = !resTable$technical_excluded & !is.na( resTable$padj ) & is.finite( resTable$stat )
        fwrite( resTable, file = paste0( tablePath, "DESeq2_", cohorts[ cohort ], "_PD_vs_HC_", population, ".txt.gz" ), sep = "\t" )

        designTab = rbind( designTab, data.frame( entry = id, cohort = cohort, stratum = population, model_formula = "~ Age_z + Diagnosis",
                                                  n_HC = sum( d$Diagnosis == "HC" ), n_PD = sum( d$Diagnosis == "PD" ), residual_df = nrow( d ) - 3,
                                                  genes_fitted = nrow( resTable ), TF_input_genes = sum( resTable$TF_input ) ) )

        # VIPER on the Wald statistics
        geneStats = data.frame( gene = toupper( trimws( resTable$gene[ resTable$TF_input ] ) ), stat = resTable$stat[ resTable$TF_input ] )
        mat = matrix( geneStats$stat, ncol = 1, dimnames = list( geneStats$gene, id ) )
        viperRes = as.data.frame( run_viper( mat, net, minsize = 5, pleiotropy = TRUE ) )

        tab = data.frame( source = names( targetsPerTF ),
                          network_targets = sapply( targetsPerTF, length ),
                          mapped_targets = sapply( targetsPerTF, function( g ) sum( g %in% geneStats$gene ) ) )
        cand = match( tab$source, viperRes$source )
        tab$NES = viperRes$score[ cand ]
        tab$PValue = viperRes$p_value[ cand ]
        tab$test_status = ifelse( tab$mapped_targets >= 5, "TESTED", "NOT_TESTED_MAPPED_TARGETS_LT5" )
        tab$FDR_BH_within_population = NA
        tab$FDR_BH_within_population[ tab$test_status == "TESTED" ] = p.adjust( tab$PValue[ tab$test_status == "TESTED" ], method = "BH" )
        tab$cohort = cohort
        tab$population = population
        tab$entry = id
        tfTable = rbind( tfTable, tab )
    }
}

fwrite( designTab, file = paste0( tablePath, "ExternalPDCohorts_Donors_And_Genes.txt" ), sep = "\t" )


#
### Supplementary Table S1
#

tested = tfTable$test_status == "TESTED"
tableS1 = data.frame( dataset = cohorts[ tfTable$cohort ], contrast = "PD vs HC", population = tfTable$population, TF = tfTable$source,
                      VIPER_NES = tfTable$NES, p_value = tfTable$PValue, FDR_BH_within_population = tfTable$FDR_BH_within_population,
                      `significant_FDR_lt_0.05` = ifelse( tested & tfTable$FDR_BH_within_population < 0.05, "TRUE", "FALSE" ),
                      test_status = tfTable$test_status, regulon_targets = tfTable$network_targets, mapped_targets = tfTable$mapped_targets,
                      model = "~ Age_z + Diagnosis", source_entry = tfTable$entry, check.names = FALSE )
tableS1$`significant_FDR_lt_0.05`[ is.na( tableS1$`significant_FDR_lt_0.05` ) ] = "FALSE"

fwrite( tableS1, file = paste0( tablePath, "Supplementary_Table_S1_TF_activity.txt" ), sep = "\t", quote = FALSE, na = "" )

# overview per data set
overview = c()
for( ds in cohorts )
{
    x = tableS1[ tableS1$dataset == ds, ]
    pops = unique( x$population )
    overview = rbind( overview, data.frame( dataset = ds, contrast = "PD vs HC", populations = length( pops ), population_labels = paste( pops, collapse = ";" ),
                                            TFs_per_population = length( targetsPerTF ), TF_population_rows = nrow( x ), tested_rows = sum( x$test_status == "TESTED" ),
                                            not_tested_rows = sum( x$test_status != "TESTED" ), significant_rows_FDR_lt_0.05 = sum( x$`significant_FDR_lt_0.05` == "TRUE" ),
                                            distinct_significant_TFs = length( unique( x$TF[ x$`significant_FDR_lt_0.05` == "TRUE" ] ) ), check.names = FALSE ) )
}

fwrite( overview, file = paste0( tablePath, "Supplementary_Table_S1_dataset_summary.txt" ), sep = "\t", quote = FALSE )
