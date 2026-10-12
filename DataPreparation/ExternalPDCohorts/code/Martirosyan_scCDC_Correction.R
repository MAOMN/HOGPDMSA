library(scCDC)
library(Seurat)
library(Matrix)
library(data.table)
library(readxl)

# change to your project directory
home_directory <- "path/to/project_directory/"

# folder with the downloaded raw data (not part of this repository, see Data/RawDataInfo.txt)
rawDataPath <- "path/to/raw_data/"

# data folder
dataPath = paste0( home_directory, "Data/" )

# output folder
outputPath = paste0( home_directory, "DataPreparation/ExternalPDCohorts/" )
codePath = paste0( outputPath, "code/" )
tablePath = paste0( outputPath, "tables/" )

set.seed( 42 )


# snRNA-seq data set of Martirosyan et al. (substantia nigra, PD and control donors), files of GEO series GSE243639
geoPath = paste0( rawDataPath, "ExternalPDCohorts/Martirosyan_GSE243639/" )

# read donor table (the table starts after 5 lines of text)
# columns: N, Brain Bank ID, Sample ID, Clinical diagnosis, Age, Sex, PMI hours, RIN measure, ...
clinical = as.data.frame( fread( paste0( geoPath, "GSE243639_Clinical_data.csv.gz" ), skip = 5, sep = ";", header = TRUE ) )
donors = data.frame( sample_id = clinical[ , "Sample ID" ],
                     diagnosis = ifelse( clinical[ , "Clinical diagnosis" ] == "Parkinson's", "PD", "HC" ),
                     sex = ifelse( tolower( clinical$Sex ) == "female", "F", "M" ),
                     age = clinical$Age )
donors = donors[ order( donors$sample_id ), ]
table( donors$diagnosis )

# read cell types of the authors
# columns: CELL_ID, UMAP_1, UMAP_2, CLUSTER, IDENT
cells = as.data.frame( read_excel( paste0( geoPath, "GSE243639_UMAP_coordinates.xlsx" ), sheet = "HighLevelCellTypes", col_types = "text" ) )
cells = data.frame( cell_id = cells$CELL_ID, sample_id = sub( "_.*$", "", cells$CELL_ID ), cell_type = cells$IDENT )
cells$diagnosis = donors$diagnosis[ match( cells$sample_id, donors$sample_id ) ]
table( cells$cell_type )

# count table (genes x cells) of all donors, the cell names end with .1 instead of -1
countFile = paste0( geoPath, "GSE243639_Filtered_count_table.csv.gz" )
countCells = colnames( fread( countFile, nrows = 0 ) )[ -1 ]
countCells = sub( "[.]1$", "-1", countCells )
length( countCells )
all( countCells %in% cells$cell_id )
cells = cells[ match( countCells, cells$cell_id ), ]


#
### Correction of ambient RNA contamination with scCDC, separately for each donor
#

corrected = list()
geneTab = list()
overview = c()
for( donor in donors$sample_id )
{
    # the count table is too large to be read at once, read the cells of one donor
    cols = which( cells$sample_id == donor )
    x <- fread( countFile, select = c( 1, cols + 1 ), showProgress = FALSE )
    counts <- Matrix( as.matrix( x[ , -1 ] ), sparse = TRUE )
    dimnames( counts ) <- list( x[[ 1 ]], cells$cell_id[ cols ] )
    rm( x )

    # all cell types of the donor are used, the cell types of the authors are the identities
    cellInfo = cells[ cols, ]
    rownames( cellInfo ) = cellInfo$cell_id
    obj <- CreateSeuratObject( counts = counts, meta.data = cellInfo, project = paste0( "Martirosyan_", donor ) )
    Idents( obj ) <- factor( obj$cell_type )
    print( paste( donor, ncol( obj ), "cells" ) )

    #contaminating genes (GCGs): genes that are found in all cell types although they should be restricted to some
    detection <- ContaminationDetection( seuratobject = obj, restriction_factor = 0.5, sample_name = paste0( "Martirosyan__", donor ), min.cell = 100, percent.cutoff = 0.2,
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
                                            raw_total_counts = sum( counts@x ), corrected_total_counts = sum( correctedCounts@x ) ) )

    # only the OPCs and oligodendrocytes are analyzed later
    corrected[[ donor ]] = correctedCounts[ , cellInfo$cell_id[ cellInfo$cell_type %in% c( "OPC", "Oligo" ) ] ]

    rm( obj, result, counts, correctedCounts )
    gc()
}
overview

fwrite( rbindlist( geneTab, fill = TRUE ), file = paste0( tablePath, "Martirosyan_scCDC_ContaminatingGenes.txt" ), sep = "\t" )
fwrite( overview, file = paste0( tablePath, "Martirosyan_scCDC_Overview.txt" ), sep = "\t" )


#
### Corrected counts of the OPCs and oligodendrocytes of all donors
#

counts <- do.call( cbind, corrected )
dim( counts )

olCells = cells[ match( colnames( counts ), cells$cell_id ), c( "cell_id", "sample_id", "diagnosis", "cell_type" ) ]
olCells$labels_major = ifelse( olCells$cell_type == "OPC", "OPCs", "Oligodendrocytes" )
rownames( olCells ) = NULL
table( olCells$labels_major, olCells$diagnosis )

# the count table of the authors has gene symbols
genes = data.frame( feature = rownames( counts ), gene_symbol = rownames( counts ) )

# list with: counts (genes x cells), cells (cell_id, sample_id, diagnosis, cell_type, labels_major), genes (feature, gene_symbol), donors (sample_id, diagnosis, sex, age)
dir.create( paste0( dataPath, "ExternalPDCohorts/OPC_OL_Counts/" ), recursive = TRUE, showWarnings = FALSE )
saveRDS( list( counts = counts, cells = olCells, genes = genes, donors = donors ), file = paste0( dataPath, "ExternalPDCohorts/OPC_OL_Counts/Martirosyan.rds" ) )
