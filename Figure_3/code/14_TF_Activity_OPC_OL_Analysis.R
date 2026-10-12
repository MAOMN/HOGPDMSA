library(decoupleR)
library(viper)
library(DESeq2)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
hogPath = paste0( home_directory, "Figure_2/" )
outputPath = paste0( home_directory, "Figure_3/" )
tablePath = paste0( outputPath, "tables/TF_Activity/" )
objectPath = paste0( outputPath, "objects/" )


# same analysis as for the HOG cells (Figure 2F), now for each cell state of the patients: MSA vs HC and PD vs HC

# signed brainSCOPE oligodendrocyte network (columns: source, target, mor, likelihood), made by 05_TF_Activity_Analysis.R of Figure 2,
# that script needs the brainSCOPE network file that has to be downloaded first
net = as.data.frame( fread( paste0( hogPath, "tables/TF_Activity/brainSCOPE_Oli_Network_SignedLinks.txt.gz" ) ) )
targetsPerTF = split( net$target, net$source )

# technical genes (ribosomal and mitochondrial genes, see Data/DataFileInfo.txt; columns: gene_symbol, ensembl_id)
technicalGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/TechnicalGenes.txt" ) ) )

# DESeq2 results (all genes) of each cell state and comparison (08_Pseudobulk_DESeq2_Analysis.R)
deseqResults <- readRDS( paste0( objectPath, "Pseudobulk_DESeq2_Results.rds" ) )

states = c( "OPC", paste0( "OL", 1:6 ) )


#
### Transcription factor activities with VIPER
#

tfTable = c()
geneNumbers = c()
for( disease in c( "MSA", "PD" ) )
{
    for( state in states )
    {
        res = as.data.frame( deseqResults[[ paste0( disease, "_vs_HC__", state ) ]] )

        # Wald statistic of all genes with adjusted p-value, without technical genes
        keep = !toupper( rownames( res ) ) %in% technicalGenes$gene_symbol & !is.na( res$padj ) & is.finite( res$stat )
        geneStats = data.frame( gene = toupper( trimws( rownames( res )[ keep ] ) ), stat = res$stat[ keep ] )
        geneNumbers = rbind( geneNumbers, data.frame( context = paste0( disease, "_vs_HC__", state ), genes = nrow( geneStats ), model = "~Age_z+Diagnosis" ) )

        mat = matrix( geneStats$stat, ncol = 1, dimnames = list( geneStats$gene, paste0( disease, "_vs_HC__", state ) ) )
        viperRes = as.data.frame( run_viper( mat, net, minsize = 5, pleiotropy = TRUE ) )

        tab = data.frame( source = names( targetsPerTF ),
                          network_targets = sapply( targetsPerTF, length ),
                          mapped_targets = sapply( targetsPerTF, function( g ) sum( g %in% geneStats$gene ) ) )
        cand = match( tab$source, viperRes$source )
        tab$score = viperRes$score[ cand ]
        tab$PValue = viperRes$p_value[ cand ]
        tab$status = ifelse( tab$mapped_targets >= 5, "TESTED", "NOT_TESTED_MAPPED_TARGETS_LT5" )
        tab$FDR = NA
        tab$FDR[ tab$status == "TESTED" ] = p.adjust( tab$PValue[ tab$status == "TESTED" ], method = "BH" )
        tab$disease = disease
        tab$stratum = state

        tfTable = rbind( tfTable, tab )
    }
}
tfTable <- tfTable[ order( tfTable$stratum, tfTable$disease, tfTable$source ), ]
rownames( tfTable ) <- NULL

fwrite( tfTable, file = paste0( tablePath, "TF_Activity_OPC_OL_AllTFs.txt" ), sep = "\t" )
fwrite( geneNumbers, file = paste0( tablePath, "TF_Activity_OPC_OL_NumberOfGenes.txt" ), sep = "\t" )

# overview
overview = c()
for( disease in c( "MSA", "PD" ) )
{
    d = tfTable[ tfTable$disease == disease, ]
    overview = rbind( overview, data.frame( disease = disease, contexts = length( unique( d$stratum ) ), tested = sum( d$status == "TESTED" ),
                                            significant_tests = sum( d$FDR < 0.05, na.rm = TRUE ),
                                            significant_regulators = length( unique( d$source[ !is.na( d$FDR ) & d$FDR < 0.05 ] ) ) ) )
}
fwrite( overview, file = paste0( tablePath, "TF_Activity_OPC_OL_Overview.txt" ), sep = "\t" )


#
### What do the transcription factors that are significant in the HOG cells do in the patients?
#

# significant transcription factors of the HOG cells (05_TF_Activity_Analysis.R of Figure 2)
# columns: source, network_targets, mapped_targets, score, PValue, status, FDR
hogTF = as.data.frame( fread( paste0( hogPath, "tables/TF_Activity/TF_Activity_HOG_SignificantTFs.txt" ) ) )

hogInPatients = tfTable[ tfTable$source %in% hogTF$source, ]
cand = match( hogInPatients$source, hogTF$source )
hogInPatients$HOG_score = hogTF$score[ cand ]
hogInPatients$HOG_FDR = hogTF$FDR[ cand ]
# significant in the patients and same direction as in the HOG cells
hogInPatients$same_direction_significant = !is.na( hogInPatients$FDR ) & hogInPatients$FDR < 0.05 & sign( hogInPatients$score ) == sign( hogInPatients$HOG_score )
hogInPatients <- hogInPatients[ order( hogInPatients$source ), ]

fwrite( hogInPatients, file = paste0( tablePath, "TF_Activity_OPC_OL_HOG_SignificantTFs.txt" ), sep = "\t" )


#
### Ranking of the transcription factors for the plot (Figure 3G): the 10 that are most consistently significant over the cell states in PD
#

# the table is read again from the saved file
tfTable = as.data.frame( fread( paste0( tablePath, "TF_Activity_OPC_OL_AllTFs.txt" ) ) )

# ranking: sum of the FDR values of the 7 cell states (PD vs HC), TFs that could not be tested in all cell states are not ranked
pd = tfTable[ tfTable$disease == "PD", ]
pdStates = unique( pd$stratum )
fdrMat = sapply( pdStates, function( s ) pd$FDR[ pd$stratum == s ][ match( names( targetsPerTF ), pd$source[ pd$stratum == s ] ) ] )
rownames( fdrMat ) = names( targetsPerTF )
colnames( fdrMat ) = paste0( "FDR_", pdStates )

ranking = data.frame( source = rownames( fdrMat ), FDR_sum = rowSums( fdrMat ), eligible = rowSums( is.na( fdrMat ) ) == 0,
                      missing_groups = apply( fdrMat, 1, function( x ) paste( pdStates[ is.na( x ) ], collapse = ";" ) ), fdrMat )
ranking <- ranking[ order( !ranking$eligible, ranking$FDR_sum, ranking$source ), ]
ranking$display_rank = NA
ranking$display_rank[ ranking$eligible ] = 1:sum( ranking$eligible )

fwrite( ranking, file = paste0( tablePath, "TF_Activity_OPC_OL_Ranking_PD_FDRSum.txt" ), sep = "\t" )
