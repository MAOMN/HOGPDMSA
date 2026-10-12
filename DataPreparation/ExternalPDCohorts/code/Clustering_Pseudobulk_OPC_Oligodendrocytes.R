library(Seurat)
library(harmony)
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


# OPC and oligodendrocyte populations in three published snRNA-seq data sets (PD and control donors)
# the cells of each data set are clustered with the steps that are used for our own data (03_OPC_OL_Clustering_Analysis.R of Figure 3),
# then the counts of each cell population are summed up per donor (pseudobulk)

cohorts = c( martirosyan = "Martirosyan", smajic = "Smajic", lee = "Lee" )

# read technical genes (ribosomal and mitochondrial genes, see Data/DataFileInfo.txt)
# columns: gene_symbol, ensembl_id (genes without symbol are named by their Ensembl id in the Smajic data set)
technicalGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/TechnicalGenes.txt" ) ) )
technicalNames = c( technicalGenes$gene_symbol, technicalGenes$ensembl_id[ technicalGenes$ensembl_id != "" ] )

# two marker genes per cell type, a cell is removed if both markers of one cell type are highly expressed
markerPairs = data.frame( type = c( "neuron", "astrocyte", "microglia", "endothelial", "pericyte" ),
                          marker_1 = c( "SNAP25", "AQP4", "P2RY12", "CLDN5", "PDGFRB" ),
                          marker_2 = c( "SYT1", "GFAP", "CX3CR1", "PECAM1", "RGS5" ) )

for( cohort in names( cohorts ) )
{
    print( cohorts[ cohort ] )

    # corrected counts of the OPCs and oligodendrocytes (cell types of the authors)
    # list with: counts (genes x cells), cells (cell_id, sample_id, diagnosis, cell_type, labels_major), genes (feature, gene_symbol), donors (sample_id, diagnosis, sex, age)
    x <- readRDS( paste0( dataPath, "ExternalPDCohorts/OPC_OL_Counts/", cohorts[ cohort ], ".rds" ) )
    cellInfo = x$cells
    rownames( cellInfo ) = cellInfo$cell_id

    obj <- CreateSeuratObject( counts = x$counts, meta.data = cellInfo, project = cohorts[ cohort ] )
    obj <- NormalizeData( obj )


    #
    ### Remove cells that also express marker genes of other cell types
    #

    # not all marker genes are in all data sets, a cell type is only checked if both of its markers are present
    markerGenes = unique( c( markerPairs$marker_1, markerPairs$marker_2 ) )
    markerGenes = markerGenes[ markerGenes %in% rownames( obj ) ]

    expr <- GetAssayData( obj, layer = "data" )[ markerGenes, ]
    rawCounts <- GetAssayData( obj, layer = "counts" )[ markerGenes, ]

    # highly expressed: above the 95% quantile of the cells that express the marker (and at least 0.25)
    isHigh <- sapply( markerGenes, function( g ) {
      v <- as.numeric( expr[ g, ] )
      expressed <- as.numeric( rawCounts[ g, ] ) > 0
      threshold <- max( quantile( v[ expressed ], 0.95, names = FALSE ), 0.25 )
      v >= threshold & expressed } )

    contaminated <- rep( FALSE, ncol( obj ) )
    for( i in 1:nrow( markerPairs ) )
    {
        if( !all( c( markerPairs$marker_1[ i ], markerPairs$marker_2[ i ] ) %in% markerGenes ) ) next
        both = isHigh[ , markerPairs$marker_1[ i ] ] & isHigh[ , markerPairs$marker_2[ i ] ]
        print( paste( markerPairs$type[ i ], sum( both ) ) )
        contaminated = contaminated | both
    }
    print( table( contaminated ) )

    obj <- subset( obj, cells = colnames( obj )[ !contaminated ] )


    #
    ### Clustering
    #

    obj <- FindVariableFeatures( obj, nfeatures = 3000 )
    varGenes <- VariableFeatures( obj )

    # technical genes and hemoglobin genes are not used for clustering, also not the marker genes of the other cell types
    isTechnical <- toupper( varGenes ) %in% technicalNames | grepl( "^HB[^(P)]", varGenes, ignore.case = TRUE )
    clusterGenes <- setdiff( varGenes[ !isTechnical ], markerGenes )
    print( length( clusterGenes ) )

    obj <- ScaleData( obj, features = clusterGenes )
    obj <- RunPCA( obj, features = clusterGenes )

    # donors are integrated with harmony
    set.seed( 42 )
    obj <- RunHarmony( obj, group.by.vars = "sample_id", dims.use = 1:30, nclust = 30, project.dim = FALSE )

    obj <- FindNeighbors( obj, reduction = "harmony", dims = 1:30 )
    obj <- FindClusters( obj, resolution = 0.3, random.seed = 42 )

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
        clusterTab = rbind( clusterTab, data.frame( natural_cluster = cl,
                                                    n_cells = nrow( m ),
                                                    n_opc_author = sum( m$labels_major == "OPCs" ),
                                                    opc_fraction = mean( m$labels_major == "OPCs" ),
                                                    n_ol_author = sum( m$labels_major == "Oligodendrocytes" ),
                                                    n_donors = length( unique( m$sample_id ) ),
                                                    n_donors_with_min_cells = sum( cellsPerDonor >= 100 ) ) )
    }
    clusterTab$final_retained <- clusterTab$n_cells >= 200 & clusterTab$n_donors_with_min_cells >= 3

    # names of the clusters: the cluster with mostly OPCs (author annotation) is OPC,
    # the oligodendrocyte clusters are numbered by size (OL1 = largest)
    clusterTab$final_cluster_label <- ifelse( !clusterTab$final_retained, "not_kept", ifelse( clusterTab$opc_fraction >= 0.5, "OPC", "OL" ) )
    olClusters <- which( clusterTab$final_cluster_label == "OL" )
    olClusters <- olClusters[ order( -clusterTab$n_cells[ olClusters ] ) ]
    clusterTab$final_cluster_label[ olClusters ] <- paste0( "OL", 1:length( olClusters ) )
    print( clusterTab )

    fwrite( clusterTab, file = paste0( tablePath, cohorts[ cohort ], "_OPC_OL_ClusterOverview.txt" ), sep = "\t" )

    meta$final_cluster_label = clusterTab$final_cluster_label[ match( meta$natural_cluster, clusterTab$natural_cluster ) ]
    meta = meta[ meta$final_cluster_label != "not_kept", ]


    #
    ### Pseudobulk counts per donor for each cell population
    #

    # gene names of the pseudobulk tables: gene symbols in upper case, genes with the same symbol are summed up
    # genes without symbol and genes with an underscore in the symbol (e.g. Y_RNA, several copies) keep their name of the count matrix
    geneNames = toupper( trimws( x$genes$gene_symbol ) )
    keepName = geneNames == "" | grepl( "_", geneNames )
    geneNames[ keepName ] = x$genes$feature[ keepName ]

    dir.create( paste0( dataPath, "ExternalPDCohorts/Pseudobulk/", cohort ), recursive = TRUE, showWarnings = FALSE )
    donorTab = c()
    for( population in sort( unique( meta$final_cluster_label ) ) )
    {
        popCells = meta[ meta$final_cluster_label == population, ]
        popDonors = sort( unique( popCells$sample_id ), method = "radix" )

        pb = sapply( popDonors, function( donor ) rowSums( x$counts[ , popCells$cell_id[ popCells$sample_id == donor ], drop = FALSE ] ) )
        pb = rowsum( pb, group = geneNames, reorder = FALSE )

        coldata = x$donors[ match( popDonors, x$donors$sample_id ), c( "sample_id", "diagnosis", "sex", "age" ) ]
        coldata$n_cells = as.integer( table( popCells$sample_id )[ popDonors ] )
        rownames( coldata ) = NULL
        donorTab = rbind( donorTab, data.frame( population = population, coldata ) )

        # list with: counts (genes x donors), coldata (sample_id, diagnosis, sex, age, n_cells)
        saveRDS( list( counts = pb, coldata = coldata ), file = paste0( dataPath, "ExternalPDCohorts/Pseudobulk/", cohort, "/", population, ".rds" ) )
    }

    # number of cells of each donor in each cell population
    fwrite( donorTab, file = paste0( tablePath, cohorts[ cohort ], "_OPC_OL_CellsPerDonorAndPopulation.txt" ), sep = "\t" )

    rm( obj, x )
    gc()
}
