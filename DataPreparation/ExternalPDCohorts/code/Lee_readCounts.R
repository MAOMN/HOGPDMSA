library(Matrix)
library(data.table)
library(rhdf5)

# reads the count matrices of the Lee data set that were made by Lee_Counts_QC_Scrublet.py (one file per donor, only cells that passed the quality control)
# needs: dataPath
# gives: counts (genes x cells, all donors), cells (cell_id, sample_id, diagnosis), genes (feature, gene_symbol, ensembl_id), donors (gsm_id, sample_id, diagnosis, sex, age)

donors = as.data.frame( fread( paste0( dataPath, "ExternalPDCohorts/Lee_Donors.txt" ) ) )

counts = list()
for( donor in donors$sample_id )
{
    h5File = paste0( dataPath, "ExternalPDCohorts/Lee_Counts/", donor, "_Counts.h5" )
    m <- sparseMatrix( i = as.integer( h5read( h5File, "matrix/indices" ) ) + 1L, p = as.integer( h5read( h5File, "matrix/indptr" ) ),
                       x = as.numeric( h5read( h5File, "matrix/data" ) ), dims = as.integer( h5read( h5File, "matrix/shape" ) ) )
    colnames( m ) <- as.character( h5read( h5File, "cells/cell_id" ) )
    counts[[ donor ]] = m
}
counts <- do.call( cbind, counts )

# the genes are the same in all files
genes = data.frame( gene_symbol = as.character( h5read( h5File, "features/name" ) ), ensembl_id = as.character( h5read( h5File, "features/id" ) ) )

# names of the genes in the count matrix: gene symbols, made unique (some symbols belong to several Ensembl ids),
# Seurat does not allow underscores in gene names, they are replaced by dashes
genes$feature = gsub( "_", "-", make.unique( genes$gene_symbol ) )
genes = genes[ , c( "feature", "gene_symbol", "ensembl_id" ) ]
rownames( counts ) <- genes$feature

cells = data.frame( cell_id = colnames( counts ), sample_id = sub( "__.*$", "", colnames( counts ) ) )
cells$diagnosis = donors$diagnosis[ match( cells$sample_id, donors$sample_id ) ]
rownames( cells ) = cells$cell_id
