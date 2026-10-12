library(Seurat)
library(Matrix)
library(reticulate)

# reads the snRNA-seq data set of the authors (all cell types, HC, MSA and PD donors) and
# provides the raw counts of single donors as Seurat objects
# needs: rawDataPath (folder with the downloaded raw data), pythonPath (Python with the package anndata)

# the raw counts of the h5ad file have no gene names, the genes are in the same order as in the normalized matrix
use_python( pythonPath, required = TRUE )
anndata <- import( "anndata", convert = FALSE )
adata <- anndata$read_h5ad( paste0( rawDataPath, "snRNAseq_HC_MSA_PD/adata.h5ad" ) )
raw <- adata$raw$to_adata()
genes = as.character( py_to_r( adata$var_names$tolist() ) )

# cell annotation of the authors: donor, diagnosis and cell type
obs <- py_to_r( adata$obs )
obs$cell_id <- rownames( obs )
diagnosisNames = c( "normal" = "HC", "multiple system atrophy" = "MSA", "Parkinson disease" = "PD" )
obs <- obs[ obs$disease__ontology_label %in% names( diagnosisNames ), ]
obs$diagnosis <- unname( diagnosisNames[ as.character( obs$disease__ontology_label ) ] )
obs$donor <- as.character( obs$donor_id )


####################################################################################################
#
# Seurat object with the raw counts of one donor, cell types of the authors as identities
#
# donor - id of the donor
#
####################################################################################################
getDonorObject <- function( donor )
{
    cells = obs$cell_id[ obs$donor == donor ]

    #counts from Python are cells x genes
    x <- raw[ cells, ]$X$tocoo()
    counts <- sparseMatrix( i = as.integer( py_to_r( x$row ) ) + 1L, j = as.integer( py_to_r( x$col ) ) + 1L, x = as.numeric( py_to_r( x$data ) ),
                            dims = c( length( cells ), length( genes ) ) )
    counts <- t( counts )
    rownames( counts ) <- genes
    colnames( counts ) <- cells

    cellInfo = data.frame( cell_id = cells, sample_id = donor, diagnosis = obs[ cells, "diagnosis" ],
                           labels_major = as.character( obs[ cells, "labels_major" ] ), labels_minor = as.character( obs[ cells, "labels_minor" ] ), row.names = cells )
    obj <- CreateSeuratObject( counts = counts, meta.data = cellInfo, project = paste0( "NatComm_", donor ) )

    #cells without cell type are not used
    obj <- subset( obj, cells = colnames( obj )[ !is.na( obj$labels_major ) ] )
    Idents( obj ) <- factor( as.character( obj$labels_major ) )
    obj[ , !is.na( Idents( obj ) ) ]
}
