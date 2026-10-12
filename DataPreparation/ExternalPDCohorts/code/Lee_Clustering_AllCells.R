library(Seurat)
library(Matrix)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/" )

# output folder
outputPath = paste0( home_directory, "DataPreparation/ExternalPDCohorts/" )
codePath = paste0( outputPath, "code/" )
tablePath = paste0( outputPath, "tables/" )

set.seed( 148434 )


# snRNA-seq data set of Lee et al. (GEO series GSE148434): there is no cell type annotation of the authors in GEO,
# so all cells are clustered here following the description of the authors and the clusters are assigned to cell types by marker genes (see Lee_scCDC_Correction.R)

# read count matrix of all donors
# gives: counts (genes x cells), cells (cell_id, sample_id, diagnosis), genes (feature, gene_symbol, ensembl_id), donors
source( paste0( codePath, "Lee_readCounts.R" ) )
dim( counts )
table( cells$sample_id )


#
### Integration of the donors and clustering
#

# mitochondrial genes are not used
obj <- CreateSeuratObject( counts = counts[ !grepl( "^MT-", rownames( counts ) ), ], meta.data = cells )
objs <- SplitObject( obj, split.by = "sample_id" )
rm( obj )
gc()

# 5000 variable genes per donor, 2000 genes for the integration
for( i in 1:length( objs ) )
{
    objs[[ i ]] <- NormalizeData( objs[[ i ]] )
    objs[[ i ]] <- FindVariableFeatures( objs[[ i ]], nfeatures = 5000 )
}
integrationGenes <- SelectIntegrationFeatures( objs, nfeatures = 2000 )
for( i in 1:length( objs ) )
{
    objs[[ i ]] <- ScaleData( objs[[ i ]], features = integrationGenes )
    objs[[ i ]] <- RunPCA( objs[[ i ]], features = integrationGenes, npcs = 45, seed.use = 148434 )
}

# reciprocal PCA, the two donors with the most cells are the references
references <- order( sapply( objs, ncol ), decreasing = TRUE )[ 1:2 ]
names( objs )[ references ]
anchors <- FindIntegrationAnchors( object.list = objs, anchor.features = integrationGenes, dims = 1:45, reduction = "rpca", reference = references )
obj <- IntegrateData( anchors, dims = 1:45 )
rm( anchors, objs )
gc()

# 45 principal components and resolution 1.5 as in the publication
DefaultAssay( obj ) <- "integrated"
obj <- ScaleData( obj )
obj <- RunPCA( obj, npcs = 45, seed.use = 148434 )
obj <- FindNeighbors( obj, reduction = "pca", dims = 1:45 )
obj <- FindClusters( obj, resolution = 1.5, random.seed = 148434 )

clusterTab = data.frame( cell_id = colnames( obj ), sample_id = obj$sample_id, cluster = as.character( obj$seurat_clusters ) )
table( clusterTab$cluster )

# cluster of each cell
fwrite( clusterTab, file = paste0( tablePath, "Lee_AllCells_Clusters.txt.gz" ), sep = "\t" )


#
### Expression of marker genes in the clusters
#

# marker genes of the oligodendrocyte lineage (OPC, committed precursors, newly formed, myelin forming and mature oligodendrocytes) and of the other cell types
markers = list( OPC = c( "PDGFRA", "CSPG4" ), COP = c( "NEU4", "SOX6", "BMP4", "GPR17", "TNS3", "FYN" ), NFOL = c( "CASR", "ITPR2" ),
                MFOL = c( "MAL", "MOG", "PLP1", "OPALIN", "SERINC5", "CTPS1", "CTPS2" ), MOL = c( "KLK6", "APOD", "TF", "PMP22", "GRM3", "JPH4" ),
                Neuron = c( "RBFOX3", "SYT1", "SNAP25", "STMN2" ), Astrocyte = c( "AQP4", "ALDH1L1", "SLC1A3", "GFAP" ), Microglia = c( "P2RY12", "CSF1R", "C1QA", "AIF1" ),
                Endothelial = c( "CLDN5", "FLT1", "VWF", "PECAM1" ), T_cell = c( "CD3D", "CD3E", "TRBC1", "TRBC2" ), Pericyte = c( "RGS5", "PDGFRB", "ACTA2", "TAGLN" ) )

# normalized expression (counts per 10000, log) from the raw counts
counts <- counts[ , clusterTab$cell_id ]
totalCounts = colSums( counts )
clusters = as.character( sort( unique( as.integer( clusterTab$cluster ) ) ) )

markerTab = c()
for( type in names( markers ) )
{
    for( g in markers[[ type ]] )
    {
        rows = which( genes$gene_symbol == g )
        if( length( rows ) == 0 ) next
        v = colSums( counts[ rows, , drop = FALSE ] )
        expr = log1p( v / totalCounts * 10000 )
        markerTab = rbind( markerTab, data.frame( cluster = clusters, n_cells = as.integer( table( clusterTab$cluster )[ clusters ] ), marker_of = type, gene = g,
                                                  mean_log_expression = sapply( clusters, function( cl ) mean( expr[ clusterTab$cluster == cl ] ) ),
                                                  fraction_detected = sapply( clusters, function( cl ) mean( v[ clusterTab$cluster == cl ] > 0 ) ) ) )
    }
}
head( markerTab )

fwrite( markerTab, file = paste0( tablePath, "Lee_AllCells_MarkerExpression.txt" ), sep = "\t" )
