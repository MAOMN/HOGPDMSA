library(Seurat)
library(Matrix)
library(AUCell)
library(BiocParallel)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
hogPath = paste0( home_directory, "Figure_2/" )
outputPath = paste0( home_directory, "Figure_3/" )
tablePath = paste0( outputPath, "tables/" )
objectPath = paste0( outputPath, "objects/" )


# question: are the nuclear OXPHOS genes that respond to alpha-syn in the HOG cells also changed in the oligodendrocytes of the patients?
# gene signature: leading edge genes of the OXPHOS subunit GSEA of the HOG cells (nuclear encoded genes only, 03_GSEA_MitoCarta_Analysis.R of Figure 2)
# columns: gene, rank, stat
leadingEdge = as.data.frame( fread( paste0( hogPath, "tables/GSEA_OXPHOS/GSEA_NuclearGenes_OXPHOS_LeadingEdge.txt" ) ) )
signature = sort( unique( leadingEdge$gene ) )

# Seurat object with the OPCs and oligodendrocytes of the kept clusters (03_OPC_OL_Clustering_Analysis.R) and the cell table for the UMAP coordinates
obj <- readRDS( paste0( objectPath, "OPC_OL_Seurat_KeptClusters.rds" ) )
counts <- GetAssayData( obj, layer = "counts" )

cells = as.data.frame( fread( paste0( tablePath, "OPC_OL_Clustering/OPC_OL_CellTable.txt.gz" ) ) )
cells <- cells[ match( colnames( obj ), cells$cell_id ), ]
rm( obj )
gc()

# technical genes (ribosomal and mitochondrial genes, see Data/DataFileInfo.txt; columns: gene_symbol, ensembl_id) are not ranked
technicalGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/TechnicalGenes.txt" ) ) )
counts <- counts[ !( toupper( rownames( counts ) ) %in% technicalGenes$gene_symbol ), ]

# signature genes that are measured in the snRNA-seq data (the HOG gene symbols are upper case)
cand = match( signature, toupper( rownames( counts ) ) )
signatureTab = data.frame( requested_gene = signature, feature = rownames( counts )[ cand ], present = !is.na( cand ) )
fwrite( signatureTab, file = paste0( tablePath, "AUCell/AUCell_OXPHOS_Signature_Genes.txt" ), sep = "\t" )

signatureGenes = signatureTab$feature[ signatureTab$present ]


#
### AUCell score of each cell
#

# number of detected genes per cell (needed below to check the ranking cutoff)
nDetected = diff( counts@p )

# ranking of the genes within each cell, ties are ranked randomly
set.seed( 1 )
rankings <- AUCell_buildRankings( counts, plotStats = FALSE, splitByBlocks = TRUE, BPPARAM = SerialParam() )

# AUC within the top 5% of the ranked genes
aucMaxRank = ceiling( 0.05 * nrow( counts ) )
auc <- AUCell_calcAUC( list( OXPHOS_LeadingEdge = signatureGenes ), rankings, aucMaxRank = aucMaxRank )
rm( rankings )
gc()

scores = data.frame( cell_id = colnames( counts ), AUC = as.numeric( getAUC( auc )[ "OXPHOS_LeadingEdge", colnames( counts ) ] ),
                     cells[ , c( "UMAP_1", "UMAP_2", "sample_id", "diagnosis", "labels_major", "labels_minor", "natural_cluster", "final_cluster_label" ) ] )
scores$diagnosis <- factor( scores$diagnosis, levels = c( "HC", "MSA", "PD" ) )

fwrite( scores, file = paste0( tablePath, "AUCell/AUCell_OXPHOS_Signature_Scores.txt.gz" ), sep = "\t" )


#
### Summary of the scores
#

summarizeScores <- function( d )
{
    data.frame( n_cells = nrow( d ), n_donors = length( unique( d$sample_id ) ), mean_auc = mean( d$AUC ), median_auc = median( d$AUC ),
                q25_auc = quantile( d$AUC, 0.25 ), q75_auc = quantile( d$AUC, 0.75 ) )
}

# per diagnosis
tab = c()
for( dx in unique( scores$diagnosis ) ) tab = rbind( tab, data.frame( diagnosis = dx, summarizeScores( scores[ scores$diagnosis == dx, ] ) ) )
fwrite( tab, file = paste0( tablePath, "AUCell/AUCell_OXPHOS_Signature_Summary_Diagnosis.txt" ), sep = "\t" )

# per diagnosis and cell state
tab = c()
groups = unique( scores[ , c( "diagnosis", "final_cluster_label" ) ] )
for( i in 1:nrow( groups ) )
{
    tab = rbind( tab, data.frame( diagnosis = groups$diagnosis[ i ], final_cluster_label = groups$final_cluster_label[ i ],
                                  summarizeScores( scores[ scores$diagnosis == groups$diagnosis[ i ] & scores$final_cluster_label == groups$final_cluster_label[ i ], ] ) ) )
}
fwrite( tab, file = paste0( tablePath, "AUCell/AUCell_OXPHOS_Signature_Summary_Diagnosis_CellState.txt" ), sep = "\t" )

# distribution over all cells
tab = data.frame( quantile_label = c( "min", "q01", "q05", "q25", "median", "mean", "q75", "q95", "q99", "max" ),
                  auc_value = c( min( scores$AUC ), quantile( scores$AUC, 0.01 ), quantile( scores$AUC, 0.05 ), quantile( scores$AUC, 0.25 ), median( scores$AUC ), mean( scores$AUC ),
                                 quantile( scores$AUC, 0.75 ), quantile( scores$AUC, 0.95 ), quantile( scores$AUC, 0.99 ), max( scores$AUC ) ) )
fwrite( tab, file = paste0( tablePath, "AUCell/AUCell_OXPHOS_Signature_ScoreDistribution.txt" ), sep = "\t" )

# check of the ranking cutoff: how many cells have less detected genes than the cutoff, and does the score depend on the number of detected genes?
tab = c()
for( dx in c( "ALL", "HC", "MSA", "PD" ) )
{
    cand = if( dx == "ALL" ) 1:nrow( scores ) else which( scores$diagnosis == dx )
    n = nDetected[ cand ]
    tab = rbind( tab, data.frame( diagnosis = dx, n_cells = length( cand ), aucMaxRank = aucMaxRank, aucMaxRank_fraction = 0.05,
                                  min = min( n ), q01 = quantile( n, 0.01 ), q05 = quantile( n, 0.05 ), q10 = quantile( n, 0.1 ), q25 = quantile( n, 0.25 ), median = median( n ),
                                  q75 = quantile( n, 0.75 ), q90 = quantile( n, 0.9 ), q95 = quantile( n, 0.95 ), q99 = quantile( n, 0.99 ), max = max( n ),
                                  n_below_aucMaxRank = sum( n < aucMaxRank ), fraction_below_aucMaxRank = mean( n < aucMaxRank ),
                                  spearman_auc_vs_detected_features = cor( scores$AUC[ cand ], n, method = "spearman" ) ) )
}
fwrite( tab, file = paste0( tablePath, "AUCell/AUCell_OXPHOS_Signature_DetectedGenes_Check.txt" ), sep = "\t" )
