library(decoupleR)
library(viper)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
degPath = paste0( home_directory, "Figure_1/" )
outputPath = paste0( home_directory, "Figure_2/" )
tablePath = paste0( outputPath, "tables/TF_Activity/" )
objectPath = paste0( outputPath, "objects/" )

writeTable <- function( x, file ) fwrite( x, file = paste0( tablePath, file ), sep = "\t" )


#
### Gene regulatory network of oligodendrocytes from brainSCOPE
#

# Emani et al., Science 2024; download see Data/ReferenceFiles/TF_Networks/DOWNLOAD.txt
# columns: TF, enhancer, promoter, TG, edgeWeight, method, celltype, Correlation, Regulation
grn = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/TF_Networks/Oli_GRN.txt" ) ) )
grn$source <- toupper( trimws( grn$TF ) )
grn$target <- toupper( trimws( grn$TG ) )

# mode of regulation: 1 for activating, -1 for repressing links, links with uncertain regulation can not be used
grn$mor <- ifelse( grn$Regulation == "Activating", 1, ifelse( grn$Regulation == "Repressing", -1, NA ) )
grn <- grn[ !is.na( grn$mor ) & grn$source != grn$target, ]

# the same TF - target gene link is listed several times (different enhancers), keep it once
net <- unique( grn[ , c( "source", "target", "mor" ) ] )

# remove links that are listed as activating and as repressing
linkId = paste( net$source, net$target )
net <- net[ !( linkId %in% linkId[ duplicated( linkId ) ] ), ]

net$likelihood <- 1
net <- net[ order( net$source, net$target ), ]
rownames( net ) <- NULL
writeTable( net, "brainSCOPE_Oli_Network_SignedLinks.txt.gz" )


#
### Gene statistics of the HOG cells: Wald statistic of alpha-syn PFF vs control
#

# tested genes (columns: gene, symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, significance), the network uses gene symbols
deg = as.data.frame( fread( paste0( degPath, "tables/DifferentialExpression/HOG_DEG_WithGeneSymbols.txt" ) ) )
deg <- deg[ !is.na( deg$symbol ) & deg$symbol != "", ]
deg$symbol <- toupper( deg$symbol )

# if a symbol occurs twice the entry with the larger statistic is taken
deg <- deg[ order( -abs( deg$stat ) ), ]
deg <- deg[ !duplicated( deg$symbol ), ]

geneStats = data.frame( gene = deg$symbol, stat = deg$stat )
geneStats <- geneStats[ order( geneStats$gene ), ]
rownames( geneStats ) <- NULL


#
### Transcription factor activities with VIPER
#

# TFs need at least 5 target genes among the measured genes, shared targets are corrected (pleiotropy)
mat <- matrix( geneStats$stat, ncol = 1, dimnames = list( geneStats$gene, "HOG_PFF_vs_Control" ) )
viperRes <- as.data.frame( run_viper( mat, net, minsize = 5, pleiotropy = TRUE ) )

# number of target genes per TF in the network and among the measured genes
targetsPerTF = split( net$target, net$source )
tfTable = data.frame( source = names( targetsPerTF ),
                      network_targets = sapply( targetsPerTF, length ),
                      mapped_targets = sapply( targetsPerTF, function( g ) sum( g %in% geneStats$gene ) ) )

cand = match( tfTable$source, viperRes$source )
tfTable$score = viperRes$score[ cand ]
tfTable$PValue = viperRes$p_value[ cand ]
tfTable$status = ifelse( tfTable$mapped_targets >= 5, "TESTED", "NOT_TESTED_MAPPED_TARGETS_LT5" )

tfTable$FDR = NA
tfTable$FDR[ tfTable$status == "TESTED" ] = p.adjust( tfTable$PValue[ tfTable$status == "TESTED" ], method = "BH" )
writeTable( tfTable, "TF_Activity_HOG_AllTFs.txt" )

# significant transcription factors
sigTF <- tfTable[ !is.na( tfTable$FDR ) & tfTable$FDR < 0.05, ]
sigTF <- sigTF[ order( sigTF$FDR ), ]
writeTable( sigTF, "TF_Activity_HOG_SignificantTFs.txt" )


#
### Target genes of the significant TFs: does the change of the target gene fit to the activity of the TF?
#

targetTab = net[ net$source %in% sigTF$source, ]
targetTab$HOG_score = sigTF$score[ match( targetTab$source, sigTF$source ) ]
targetTab$Wald = geneStats$stat[ match( targetTab$target, geneStats$gene ) ]
targetTab <- targetTab[ !is.na( targetTab$Wald ), ]
targetTab$supports_HOG_direction = targetTab$mor * targetTab$Wald * sign( targetTab$HOG_score ) > 0
targetTab <- targetTab[ order( targetTab$target, targetTab$source ), c( "target", "source", "mor", "likelihood", "HOG_score", "Wald", "supports_HOG_direction" ) ]
writeTable( targetTab, "TF_Activity_HOG_SignificantTFs_TargetGenes.txt.gz" )

# data of the plot in 06_TF_Activity_Plots.R
saveRDS( sigTF, file = paste0( objectPath, "TF_Activity_SignificantTFs.rds" ) )
