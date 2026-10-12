library(tradeSeq)
library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(edgeR)
library(BiocParallel)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/" )

# results of the clustering (Figure 3)
clusterPath = paste0( home_directory, "Figure_3/" )

# output folder
outputPath = paste0( home_directory, "Figure_4/" )
tablePath = paste0( outputPath, "tables/fitGAM/" )
objectPath = paste0( outputPath, "objects/" )


# slingshot result (made by 01_Trajectories_Slingshot_Analysis.R)
# list with: slingshot_sce, lineages, pseudotime, cell_weights, cell_metadata
trajectory <- readRDS( paste0( objectPath, "OPC_OL_Slingshot.rds" ) )
cells <- trajectory$cell_metadata

# Seurat object with the OPCs and oligodendrocytes of the kept clusters (made by 03_OPC_OL_Clustering_Analysis.R of Figure 3)
obj <- readRDS( paste0( clusterPath, "objects/OPC_OL_Seurat_KeptClusters.rds" ) )
counts <- GetAssayData( obj, layer = "counts" )
counts <- counts[ , cells$cell_id ]
rm( obj )
gc()


#
### Genes for the models
#

# technical genes (ribosomal and mitochondrial genes, see Data/DataFileInfo.txt; columns: gene_symbol, ensembl_id)
technicalGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/TechnicalGenes.txt" ) ) )

# a gene has to be detected in at least 100 cells, technical genes are not modelled
geneTab = data.frame( gene = rownames( counts ), detected_cells = rowSums( counts > 0 ), total_counts = rowSums( counts ) )
geneTab$technical_gene <- toupper( geneTab$gene ) %in% technicalGenes$gene_symbol
geneTab$fitted_candidate <- !geneTab$technical_gene & geneTab$detected_cells >= 100

fwrite( geneTab, file = paste0( tablePath, "fitGAM_GenesUsed.txt.gz" ), sep = "\t" )

counts <- counts[ geneTab$gene[ geneTab$fitted_candidate ], ]
saveRDS( counts, file = paste0( objectPath, "fitGAM_Counts.rds" ) )

# offset of each cell: library size with TMM normalization factor
libSize <- colSums( counts )
normFactors <- calcNormFactors( counts, lib.size = libSize, method = "TMM" )
offset <- log( libSize * normFactors )
names( offset ) <- colnames( counts )
saveRDS( offset, file = paste0( objectPath, "fitGAM_Offset_ThreeLineages.rds" ) )

# conditions
conditions <- factor( as.character( cells$diagnosis ), levels = c( "HC", "MSA", "PD" ) )

# the models are fitted in batches of 100 genes
# gene_rank: position of the gene in the count matrix, row_in_batch: position of the gene within its batch
batchTab = data.frame( gene = rownames( counts ), gene_rank = 1:nrow( counts ) )
batchTab$batch_id = ceiling( batchTab$gene_rank / 100 )
batchTab$row_in_batch = batchTab$gene_rank - ( batchTab$batch_id - 1 ) * 100


# number of knots from evaluateK (05_evaluateK_Analysis.R)
nKnots = 6


#
### Fit of the models: three lineages, each with the three conditions HC, MSA and PD
#

# fitGAM needs several days for all genes, the models are saved in batches of 100 genes
# (the article used x86_64; on arm64, e.g. Apple silicon, the same seed gives other fits)
modelPath = paste0( objectPath, "fitGAM_HC_MSA_PD_ThreeLineages/" )
dir.create( modelPath, showWarnings = FALSE )
sds <- SlingshotDataSet( trajectory$slingshot_sce )
for( b in unique( batchTab$batch_id ) )
{
    genes <- batchTab$gene[ batchTab$batch_id == b ]
    set.seed( 202608072 )
    fit <- fitGAM( counts = counts[ genes, ], sds = sds, conditions = conditions, offset = offset, nknots = nKnots,
                   parallel = TRUE, BPPARAM = MulticoreParam( 4, RNGseed = 202608072 ), control = mgcv::gam.control( maxit = 100 ), sce = TRUE )
    saveRDS( fit, file = paste0( modelPath, "fit_sce_batch_", sprintf( "%04d", b ), ".rds" ) )
}

# did the fit converge for all genes? tradeSeq stores this for each gene in the model object
batchTab$converged = NA
for( b in unique( batchTab$batch_id ) )
{
    fit <- readRDS( paste0( modelPath, "fit_sce_batch_", sprintf( "%04d", b ), ".rds" ) )
    batchTab$converged[ match( rownames( fit ), batchTab$gene ) ] = rowData( fit )$tradeSeq$converged
}

# genes with their batch and convergence of the fit
fwrite( batchTab, file = paste0( tablePath, "fitGAM_GeneBatches_ThreeLineages.txt.gz" ), sep = "\t" )
