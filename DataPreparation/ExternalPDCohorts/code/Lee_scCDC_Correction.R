library(scCDC)
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

set.seed( 42 )


# read count matrix of all donors of the Lee data set
# gives: counts (genes x cells), cells (cell_id, sample_id, diagnosis), genes (feature, gene_symbol, ensembl_id), donors
source( paste0( codePath, "Lee_readCounts.R" ) )

# read cluster of each cell
# columns: cell_id, sample_id, cluster
clusterTab = as.data.frame( fread( paste0( tablePath, "Lee_AllCells_Clusters.txt.gz" ) ) )
cells = cells[ clusterTab$cell_id, ]
cells$cluster = clusterTab$cluster

# cell types of the clusters, assigned by the expression of the marker genes (Lee_AllCells_MarkerExpression.txt)
# clusters 23, 24 and 29 express markers of several cell types (Mixed), cluster 32 has no clear marker expression (Unresolved)
cellTypes = list( Oligodendrocytes = c( 0:5, 7, 9, 13, 17, 26:28 ), OPCs = 6, Astrocytes = c( 8, 12, 20, 21, 31, 33 ), Microglia = c( 11, 14 ), Neurons = c( 15, 18, 19, 22, 25, 34 ),
                  Endothelial = 10, Pericytes = 16, "T cells" = 30, Mixed = c( 23, 24, 29 ), Unresolved = 32 )
cells$cell_type = NA
for( type in names( cellTypes ) )
{
    cells$cell_type[ cells$cluster %in% cellTypes[[ type ]] ] = type
}
table( cells$cell_type, useNA = "ifany" )


#
### Correction of ambient RNA contamination with scCDC, separately for each donor
#

corrected = list()
geneTab = list()
overview = c()
for( donor in unique( cells$sample_id ) )
{
    # all cell types of the donor are used, the cell types are the identities (in the order in which they occur in the cells of the donor)
    cellInfo = cells[ cells$sample_id == donor, ]
    donorCounts = counts[ , cellInfo$cell_id ]
    obj <- CreateSeuratObject( counts = donorCounts, meta.data = cellInfo, project = paste0( "Lee_", donor ) )
    Idents( obj ) <- cellInfo$cell_type
    print( paste( donor, ncol( obj ), "cells" ) )

    #contaminating genes (GCGs): genes that are found in all cell types although they should be restricted to some
    detection <- ContaminationDetection( seuratobject = obj, restriction_factor = 0.5, sample_name = paste0( "Lee__", donor ), min.cell = 100, percent.cutoff = 0.2,
                                         out_path.plot = NULL, out_path.table = NULL, only.cont_genes = TRUE )

    #contamination ratio of the donor
    ratio <- ContaminationQuantification( object = obj, cont_genes = rownames( detection ), auc_thres = 0.9, min.cell = 100, top10_gcg = FALSE, layer = "counts" )

    #scCDC returns the Seurat object with the additional assay Corrected
    result <- ContaminationCorrection( object = obj, cont_genes = rownames( detection ), auc_thres = 0.9, min.cell = 50 )
    correctedCounts <- GetAssayData( result, assay = "Corrected", layer = "counts" )

    #one column per cell type of the donor (with at least 100 cells) and the mean distance
    geneTab[[ donor ]] = data.frame( diagnosis = cellInfo$diagnosis[ 1 ], sample_id = donor, gene = rownames( detection ), detection, check.names = FALSE, row.names = NULL )
    overview = rbind( overview, data.frame( diagnosis = cellInfo$diagnosis[ 1 ], sample_id = donor, n_cells = ncol( obj ), n_features = nrow( obj ),
                                            n_cell_types = length( levels( Idents( obj ) ) ), n_cell_types_100_cells = sum( table( Idents( obj ) ) >= 100 ),
                                            n_gcg = nrow( detection ), contamination_ratio = as.numeric( ratio ),
                                            raw_total_counts = sum( donorCounts@x ), corrected_total_counts = sum( correctedCounts@x ) ) )

    # only the OPCs and oligodendrocytes are analyzed later
    corrected[[ donor ]] = correctedCounts[ , cellInfo$cell_id[ cellInfo$cell_type %in% c( "OPCs", "Oligodendrocytes" ) ] ]

    rm( obj, result, donorCounts, correctedCounts )
    gc()
}
overview

fwrite( rbindlist( geneTab, fill = TRUE ), file = paste0( tablePath, "Lee_scCDC_ContaminatingGenes.txt" ), sep = "\t" )
fwrite( overview, file = paste0( tablePath, "Lee_scCDC_Overview.txt" ), sep = "\t" )


#
### Corrected counts of the OPCs and oligodendrocytes of all donors
#

counts <- do.call( cbind, corrected )
dim( counts )

olCells = cells[ colnames( counts ), c( "cell_id", "sample_id", "diagnosis", "cell_type" ) ]
olCells$labels_major = olCells$cell_type
rownames( olCells ) = NULL
table( olCells$labels_major, olCells$diagnosis )

# list with: counts (genes x cells), cells (cell_id, sample_id, diagnosis, cell_type, labels_major), genes (feature, gene_symbol, ensembl_id), donors (sample_id, diagnosis, sex, age)
dir.create( paste0( dataPath, "ExternalPDCohorts/OPC_OL_Counts/" ), recursive = TRUE, showWarnings = FALSE )
saveRDS( list( counts = counts, cells = olCells, genes = genes, donors = donors[ , c( "sample_id", "diagnosis", "sex", "age" ) ] ), file = paste0( dataPath, "ExternalPDCohorts/OPC_OL_Counts/Lee.rds" ) )
