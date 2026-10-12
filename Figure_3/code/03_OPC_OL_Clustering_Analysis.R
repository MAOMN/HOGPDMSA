library(Seurat)
library(harmony)
library(Matrix)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_3/" )
tablePath = paste0( outputPath, "tables/OPC_OL_Clustering/" )
objectPath = paste0( outputPath, "objects/" )


# count matrix of all cells (genes x cells), ambient RNA contamination was corrected with scCDC
# made by DataPreparation/snRNAseq_HC_MSA_PD/code/scCDC_Correction.R from the data of Rydbirk et al. (2026)
counts <- readRDS( paste0( dataPath, "snRNAseq_HC_MSA_PD/NatComm_CorrectedCounts_scCDC.rds" ) )

# cell annotation of the authors (columns: cell_id, sample_id, diagnosis, labels_major, labels_minor, UMAP_1, UMAP_2)
cells = as.data.frame( fread( paste0( dataPath, "snRNAseq_HC_MSA_PD/NatComm_CellAnnotation.txt.gz" ) ) )
rownames( cells ) <- cells$cell_id

# only oligodendrocytes and OPCs are analyzed, the UMAP of the authors is not needed any more
olCells <- cells[ cells$labels_major %in% c( "Oligodendrocytes", "OPCs" ), 1:5 ]

obj <- CreateSeuratObject( counts = counts[ , olCells$cell_id ], meta.data = olCells )
rm( counts )
gc()

obj$diagnosis <- factor( obj$diagnosis, levels = c( "HC", "MSA", "PD" ) )
obj <- NormalizeData( obj )


#
### Remove cells that also express marker genes of other cell types
#

# two marker genes per cell type, a cell is removed if both markers of one cell type are highly expressed
markerPairs = data.frame( type = c( "neuron", "astrocyte", "microglia", "endothelial", "pericyte" ),
                          marker_1 = c( "SNAP25", "AQP4", "P2RY12", "CLDN5", "PDGFRB" ),
                          marker_2 = c( "SYT1", "GFAP", "CX3CR1", "PECAM1", "RGS5" ) )
markerGenes = unique( c( markerPairs$marker_1, markerPairs$marker_2 ) )

expr <- GetAssayData( obj, layer = "data" )[ markerGenes, ]
rawCounts <- GetAssayData( obj, layer = "counts" )[ markerGenes, ]

# highly expressed: above the 95% quantile of the cells that express the marker (and at least 0.25)
isHigh <- sapply( markerGenes, function( g ) {
  x <- as.numeric( expr[ g, ] )
  expressed <- as.numeric( rawCounts[ g, ] ) > 0
  threshold <- max( quantile( x[ expressed ], 0.95, names = FALSE ), 0.25 )
  x >= threshold & expressed } )

contaminated <- rep( FALSE, ncol( obj ) )
for( i in 1:nrow( markerPairs ) ) contaminated = contaminated | ( isHigh[ , markerPairs$marker_1[ i ] ] & isHigh[ , markerPairs$marker_2[ i ] ] )

obj <- subset( obj, cells = colnames( obj )[ !contaminated ] )


#
### Clustering
#

obj <- FindVariableFeatures( obj, nfeatures = 3000 )
varGenes <- VariableFeatures( obj )

# technical genes (ribosomal and mitochondrial genes, see Data/DataFileInfo.txt; columns: gene_symbol, ensembl_id)
technicalGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/TechnicalGenes.txt" ) ) )

# technical genes and hemoglobin genes are not used for clustering, also not the marker genes of the other cell types
isTechnical <- toupper( varGenes ) %in% technicalGenes$gene_symbol | grepl( "^HB[^(P)]", varGenes, ignore.case = TRUE )
clusterGenes <- setdiff( varGenes[ !isTechnical ], markerGenes )

obj <- ScaleData( obj, features = clusterGenes )
obj <- RunPCA( obj, features = clusterGenes )

# donors are integrated with harmony
set.seed( 42 )
obj <- RunHarmony( obj, group.by.vars = "sample_id", dims.use = 1:30, nclust = 30, project.dim = FALSE )

obj <- FindNeighbors( obj, reduction = "harmony", dims = 1:30 )
obj <- FindClusters( obj, resolution = 0.3, random.seed = 42 )
obj <- RunUMAP( obj, reduction = "harmony", dims = 1:30 )

obj$natural_cluster <- as.character( Idents( obj ) )


#
### Which clusters are kept
#

meta <- obj@meta.data
clusters = as.character( sort( unique( as.integer( meta$natural_cluster ) ) ) )

# a cluster needs at least 200 cells and at least 3 donors that have 100 or more cells in it
clusterTab = c()
for( cl in clusters )
{
    m = meta[ meta$natural_cluster == cl, ]
    cellsPerDonor = table( m$sample_id )
    donorsOK = cellsPerDonor[ cellsPerDonor >= 100 ]
    donorsOK = donorsOK[ order( -donorsOK, names( donorsOK ), method = "radix" ) ]

    clusterTab = rbind( clusterTab, data.frame( natural_cluster = cl,
                                                n_cells = nrow( m ),
                                                n_opc_author = sum( m$labels_major == "OPCs" ),
                                                opc_fraction = mean( m$labels_major == "OPCs" ),
                                                n_ol_author = sum( m$labels_major == "Oligodendrocytes" ),
                                                n_donors = length( unique( m$sample_id ) ),
                                                n_diagnoses = length( unique( m$diagnosis ) ),
                                                donor_support_min_donors = 3,
                                                donor_support_min_cells_per_donor = 100,
                                                n_donors_with_min_cells = length( donorsOK ),
                                                supported_donor_cell_counts = paste0( names( donorsOK ), ":", donorsOK, collapse = ";" ) ) )
}
clusterTab$supported_donor_cell_counts[ clusterTab$n_donors_with_min_cells == 0 ] <- ""
clusterTab$donor_support_pass <- clusterTab$n_donors_with_min_cells >= 3
clusterTab$small_cluster <- clusterTab$n_cells < 200
clusterTab$final_retained <- !clusterTab$small_cluster & clusterTab$donor_support_pass
clusterTab$final_exclusion_reason <- ifelse( clusterTab$small_cluster, "small_cluster", ifelse( !clusterTab$donor_support_pass, "donor_unsupported", "retained" ) )

# names of the clusters: the cluster with mostly OPCs (author annotation) is OPC, the oligodendrocyte clusters are numbered by size (OL1 = largest)
clusterTab$broad_label_class <- ifelse( clusterTab$final_retained & clusterTab$opc_fraction >= 0.5, "OPC",
                                        ifelse( clusterTab$final_retained, "OL", ifelse( clusterTab$small_cluster, "Small", "DonorUnsupported" ) ) )
clusterTab$final_cluster_label <- clusterTab$broad_label_class
olClusters <- which( clusterTab$broad_label_class == "OL" )
olClusters <- olClusters[ order( -clusterTab$n_cells[ olClusters ] ) ]
clusterTab$final_cluster_label[ olClusters ] <- paste0( "OL", 1:length( olClusters ) )

fwrite( clusterTab, file = paste0( tablePath, "OPC_OL_ClusterOverview.txt" ), sep = "\t" )

# add to the cells
cand = match( meta$natural_cluster, clusterTab$natural_cluster )
for( field in c( "final_cluster_label", "final_retained", "small_cluster", "donor_support_pass", "n_donors_with_min_cells", "final_exclusion_reason" ) )
{
    obj[[ field ]] <- clusterTab[ cand, field ]
}

# number of cells of each donor in each cluster
donorTab = as.data.frame( table( sample_id = obj$sample_id, natural_cluster = obj$natural_cluster ), stringsAsFactors = FALSE )
donorTab <- donorTab[ donorTab$Freq > 0, ]
colnames( donorTab )[ 3 ] <- "n_cells"
cand = match( donorTab$natural_cluster, clusterTab$natural_cluster )
donorTab$final_cluster_label <- clusterTab$final_cluster_label[ cand ]
donorTab$final_retained <- clusterTab$final_retained[ cand ]
donorTab$diagnosis <- as.character( obj$diagnosis[ match( donorTab$sample_id, obj$sample_id ) ] )
donorTab$cluster_total <- clusterTab$n_cells[ cand ]
donorTab$donor_fraction <- donorTab$n_cells / donorTab$cluster_total
donorTab$donor_ge_min_cells <- donorTab$n_cells >= 100
donorTab <- donorTab[ order( as.integer( donorTab$natural_cluster ), -donorTab$donor_fraction, donorTab$sample_id, method = "radix" ),
                      c( "natural_cluster", "final_cluster_label", "final_retained", "sample_id", "diagnosis", "n_cells", "cluster_total", "donor_fraction", "donor_ge_min_cells" ) ]

fwrite( donorTab, file = paste0( tablePath, "OPC_OL_CellsPerDonorAndCluster.txt" ), sep = "\t" )


#
### Save cell table and Seurat object
#

# table of all cells that were clustered with UMAP coordinates and cluster
umap <- Embeddings( obj, "umap" )
cellTab = data.frame( cell_id = colnames( obj ), UMAP_1 = umap[ , 1 ], UMAP_2 = umap[ , 2 ],
                      obj@meta.data[ , c( "nCount_RNA", "nFeature_RNA", "sample_id", "diagnosis", "labels_major", "labels_minor", "natural_cluster", "final_cluster_label",
                                           "final_retained", "small_cluster", "donor_support_pass", "n_donors_with_min_cells", "final_exclusion_reason" ) ] )
fwrite( cellTab, file = paste0( tablePath, "OPC_OL_CellTable.txt.gz" ), sep = "\t" )

# harmony corrected principal components of all clustered cells
saveRDS( Embeddings( obj, "harmony" ), file = paste0( objectPath, "OPC_OL_HarmonyEmbedding.rds" ) )

# Seurat object with the cells of the kept clusters only
obj <- subset( obj, cells = colnames( obj )[ obj$final_retained ] )
saveRDS( obj, file = paste0( objectPath, "OPC_OL_Seurat_KeptClusters.rds" ) )
