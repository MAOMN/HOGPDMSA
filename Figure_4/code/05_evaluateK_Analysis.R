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
tablePath = paste0( outputPath, "tables/evaluateK/" )
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
### Genes, offset and conditions
#

# technical genes (ribosomal and mitochondrial genes, see Data/DataFileInfo.txt; columns: gene_symbol, ensembl_id)
technicalGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/TechnicalGenes.txt" ) ) )

# a gene has to be detected in at least 100 cells, technical genes are not modelled
keep <- !( toupper( rownames( counts ) ) %in% technicalGenes$gene_symbol ) & rowSums( counts > 0 ) >= 100
counts <- counts[ keep, ]

# offset of each cell: library size with TMM normalization factor
libSize <- colSums( counts )
normFactors <- calcNormFactors( counts, lib.size = libSize, method = "TMM" )
offset <- log( libSize * normFactors )
names( offset ) <- colnames( counts )

# conditions
conditions <- factor( as.character( cells$diagnosis ), levels = c( "HC", "MSA", "PD" ) )


#
### evaluateK: models with 3 to 10 knots for 100 randomly selected genes, separately for each lineage
#

# cells of each lineage: all cells that belong to the lineage, with their pseudotime in this lineage
pt <- slingPseudotime( trajectory$slingshot_sce, na = FALSE )
branch <- as.character( slingBranchID( trajectory$slingshot_sce ) )

## evaluateK fits the models of the 100 genes for each number of knots (10 to 15 minutes per lineage with 8 cores), therefore the AIC values were saved (one file per lineage)
#for( i in 1:3 )
#{
#    inLineage <- grepl( paste0( "(^|,)", i, "(,|$)" ), branch )
#    time <- matrix( pt[ inLineage, i ], ncol = 1, dimnames = list( colnames( counts )[ inLineage ], paste0( "L", i ) ) )
#    weight <- matrix( 1, nrow = sum( inLineage ), ncol = 1, dimnames = dimnames( time ) )
#
#    set.seed( 202608073 )
#    aic <- evaluateK( counts = counts[ , inLineage ], k = 3:10, nGenes = 100, pseudotime = time, cellWeights = weight, conditions = conditions[ inLineage ], offset = offset[ inLineage ],
#                      plot = FALSE, verbose = FALSE, parallel = TRUE, BPPARAM = MulticoreParam( 8, RNGseed = 202608073 ), control = mgcv::gam.control( maxit = 100 ) )
#    saveRDS( list( aic = aic ), file = paste0( dataPath, "SavedModels/evaluateK/L", i, "_k3_10.rds" ) )
#}

# saved AIC values: genes x number of knots
knotTab = c()
for( i in 1:3 )
{
    aic = readRDS( paste0( dataPath, "SavedModels/evaluateK/L", i, "_k3_10.rds" ) )$aic
    k = as.integer( sub( "k: ", "", colnames( aic ) ) )

    # check that the saved AIC values belong to the genes from above (evaluateK draws the 100 genes out of the rows of the count matrix) and to 3 to 10 knots
    set.seed( 202608073 )
    selectedGenes <- rownames( counts )[ sample( seq_len( nrow( counts ) ), 100 ) ]
    print( c( lineage = i,
              cells = sum( grepl( paste0( "(^|,)", i, "(,|$)" ), branch ) ),
              genes = identical( rownames( aic ), selectedGenes ),
              knots = identical( k, 3:10 ),
              all_models_fitted = all( is.finite( aic ) ) ) )

    meanAIC = colMeans( aic )
    aicRange = apply( aic, 1, function( x ) diff( range( x ) ) )
    bestK = k[ apply( aic, 1, which.min ) ]
    knotTab = rbind( knotTab, data.frame( lineage = paste0( "L", i ), endpoint = c( "OL5", "OL6", "OL3" )[ i ], genes = nrow( aic ), k = k,
                                          mean_AIC = meanAIC, decrease_to_previous_k = c( NA, -diff( meanAIC ) ), relative_AIC = colMeans( aic / aic[ , 1 ] ),
                                          genes_with_lowest_AIC = as.integer( table( factor( bestK[ aicRange > 2 ], levels = k ) ) ) ) )
}

# mean AIC for each lineage and number of knots
# in all three lineages the decrease of the mean AIC levels off beyond 6 knots
fwrite( knotTab, file = paste0( tablePath, "fitGAM_NumberOfKnots_AIC.txt" ), sep = "\t" )
