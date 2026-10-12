library(scCDC)
library(Seurat)
library(Matrix)
library(data.table)
library(org.Hs.eg.db)

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


# snRNA-seq data set of Smajic et al. (midbrain, PD and control donors), files of GEO series GSE157783
geoPath = paste0( rawDataPath, "ExternalPDCohorts/Smajic_GSE157783/" )

# read donor information from the series file of GEO (sample titles IPD1-IPD5 and C1-C6, in the cell table the PD donors are PD1-PD5)
soft = readLines( paste0( geoPath, "GSE157783_family.soft.gz" ) )
titles = sub( "^!Sample_title = ", "", grep( "^!Sample_title", soft, value = TRUE ) )
donors = data.frame( sample_id = sub( "^IPD", "PD", titles ),
                     diagnosis = ifelse( grepl( "^IPD", titles ), "PD", "HC" ),
                     sex = sub( ".*Sex: ", "", grep( "^!Sample_characteristics_ch1 = Sex:", soft, value = TRUE ) ),
                     age = as.numeric( sub( ".*age at death: ", "", grep( "^!Sample_characteristics_ch1 = age at death:", soft, value = TRUE ) ) ) )
donors

# read cell types of the authors
# columns: barcode, cell_ontology, patient
cells = as.data.frame( fread( cmd = paste0( "tar -xzOf '", geoPath, "GSE157783_IPDCO_hg_midbrain_cell.tar.gz'" ) ) )
cells = data.frame( cell_id = cells$barcode, sample_id = cells$patient, cell_type = cells$cell_ontology )
cells$diagnosis = donors$diagnosis[ match( cells$sample_id, donors$sample_id ) ]
table( cells$cell_type )

# read gene table, the genes are given by their Ensembl ids
# columns: gene, row, ...
geneInfo = as.data.frame( fread( cmd = paste0( "tar -xzOf '", geoPath, "GSE157783_IPDCO_hg_midbrain_genes.tar.gz'" ) ) )
genes = data.frame( ensembl_id = sub( "[.][0-9]+$", "", geneInfo$gene ) )

# gene symbols (if there are several symbols for one Ensembl id the first one in alphabetical order is used)
symbols <- select( org.Hs.eg.db, keys = genes$ensembl_id, keytype = "ENSEMBL", columns = "SYMBOL" )
symbols <- symbols[ !is.na( symbols$SYMBOL ), ]
symbols <- symbols[ order( symbols$ENSEMBL, symbols$SYMBOL, method = "radix" ), ]
genes$gene_symbol = symbols$SYMBOL[ match( genes$ensembl_id, symbols$ENSEMBL ) ]
genes$gene_symbol[ is.na( genes$gene_symbol ) ] = ""
table( genes$gene_symbol != "" )

# names of the genes in the count matrix: gene symbol, Ensembl id for the genes without symbol
genes$feature = ifelse( genes$gene_symbol != "", genes$gene_symbol, genes$ensembl_id )
genes$feature = make.unique( gsub( "_", "-", genes$feature ), sep = "-dup" )

# read count table (genes x cells, without gene names, genes in the order of the gene table)
x <- fread( cmd = paste0( "tar -xzOf '", geoPath, "GSE157783_IPDCO_hg_midbrain_UMI.tar.gz'" ), colClasses = "integer" )
dim( x )
identical( colnames( x ), cells$cell_id )

# the conversion to a sparse matrix is done in blocks of 512 genes
blocks = split( 1:nrow( x ), ceiling( 1:nrow( x ) / 512 ) )
rawCounts <- do.call( rbind, lapply( blocks, function( rows ) Matrix( as.matrix( x[ rows ] ), sparse = TRUE ) ) )
dimnames( rawCounts ) <- list( genes$feature, cells$cell_id )
rm( x )
gc()


#
### Correction of ambient RNA contamination with scCDC, separately for each donor
#

corrected = list()
geneTab = list()
overview = c()
for( donor in unique( cells$sample_id ) )
{
    # all cell types of the donor are used, the cell types of the authors are the identities
    cellInfo = cells[ cells$sample_id == donor, ]
    rownames( cellInfo ) = cellInfo$cell_id
    counts = rawCounts[ , cellInfo$cell_id ]
    obj <- CreateSeuratObject( counts = counts, meta.data = cellInfo, project = paste0( "Smajic_", donor ) )
    Idents( obj ) <- factor( obj$cell_type )
    print( paste( donor, ncol( obj ), "cells" ) )

    #contaminating genes (GCGs): genes that are found in all cell types although they should be restricted to some
    detection <- ContaminationDetection( seuratobject = obj, restriction_factor = 0.5, sample_name = paste0( "Smajic__", donor ), min.cell = 100, percent.cutoff = 0.2,
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
    corrected[[ donor ]] = correctedCounts[ , cellInfo$cell_id[ cellInfo$cell_type %in% c( "OPCs", "Oligodendrocytes" ) ] ]

    rm( obj, result, counts, correctedCounts )
    gc()
}
overview

fwrite( rbindlist( geneTab, fill = TRUE ), file = paste0( tablePath, "Smajic_scCDC_ContaminatingGenes.txt" ), sep = "\t" )
fwrite( overview, file = paste0( tablePath, "Smajic_scCDC_Overview.txt" ), sep = "\t" )


#
### Corrected counts of the OPCs and oligodendrocytes of all donors, cells in the same order as in the data set of the authors
#

counts <- do.call( cbind, corrected )
olCells = cells[ cells$cell_type %in% c( "OPCs", "Oligodendrocytes" ), c( "cell_id", "sample_id", "diagnosis", "cell_type" ) ]
counts <- counts[ , olCells$cell_id ]
dim( counts )

olCells$labels_major = olCells$cell_type
rownames( olCells ) = NULL
table( olCells$labels_major, olCells$diagnosis )

# list with: counts (genes x cells), cells (cell_id, sample_id, diagnosis, cell_type, labels_major), genes (feature, gene_symbol, ensembl_id), donors (sample_id, diagnosis, sex, age)
dir.create( paste0( dataPath, "ExternalPDCohorts/OPC_OL_Counts/" ), recursive = TRUE, showWarnings = FALSE )
saveRDS( list( counts = counts, cells = olCells, genes = genes[ , c( "feature", "gene_symbol", "ensembl_id" ) ], donors = donors ), file = paste0( dataPath, "ExternalPDCohorts/OPC_OL_Counts/Smajic.rds" ) )
