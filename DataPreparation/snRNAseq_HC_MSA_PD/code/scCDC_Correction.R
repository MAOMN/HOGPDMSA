library(scCDC)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# folder with the downloaded raw data (not part of this repository, see Data/snRNAseq_HC_MSA_PD/README.txt)
rawDataPath <- "path/to/raw_data/"

# Python with the package anndata, needed to read the h5ad file (change to your Python, see PythonPackages.txt)
pythonPath <- "path/to/python"

# data folder
dataPath = paste0( home_directory, "Data/snRNAseq_HC_MSA_PD/" )

# output folder
outputPath = paste0( home_directory, "DataPreparation/snRNAseq_HC_MSA_PD/" )
codePath = paste0( outputPath, "code/" )
tablePath = paste0( outputPath, "tables/" )

set.seed( 42 )


# read the snRNA-seq data set of the authors
# gives: obs (cell annotation with donor and diagnosis), genes, function getDonorObject
source( paste0( codePath, "readAuthorCounts.R" ) )

# read contaminating genes of each donor
# columns: diagnosis, sample_id, gene, one column per cell type, mean_distance, n_cells, n_features, n_clusters
gcgTab = as.data.frame( fread( paste0( tablePath, "scCDC_ContaminatingGenes.txt" ) ) )
table( gcgTab$sample_id )


#
### Correction of the counts of the contaminating genes with scCDC, separately for each donor
#

donors = unique( obs$donor )
length( donors )

corrected = list()
overview = c()
for( donor in donors )
{
    obj <- getDonorObject( donor )
    gcg = intersect( unique( gcgTab$gene[ gcgTab$sample_id == donor ] ), rownames( obj ) )
    print( paste( donor, ncol( obj ), "cells", length( gcg ), "contaminating genes" ) )

    result <- ContaminationCorrection( object = obj, cont_genes = gcg, auc_thres = 0.9, min.cell = 50 )

    #scCDC returns the Seurat object with the additional assay Corrected
    corrected[[ donor ]] = GetAssayData( result, assay = "Corrected", layer = "counts" )
    rawCounts = GetAssayData( result, assay = "RNA", layer = "counts" )
    overview = rbind( overview, data.frame( sample_id = donor, n_cells = ncol( rawCounts ), n_features = nrow( rawCounts ), raw_total_counts = sum( rawCounts@x ),
                                            corrected_total_counts = sum( corrected[[ donor ]]@x ), nnz = length( corrected[[ donor ]]@x ) ) )
    rm( obj, result, rawCounts )
    gc()
}
overview

fwrite( overview, file = paste0( tablePath, "scCDC_Correction_Overview.txt" ), sep = "\t" )


#
### Count matrix of all cells (genes x cells), cells in the same order as in the data set of the authors
#

counts <- do.call( cbind, corrected )
counts <- counts[ , obs$cell_id ]
dim( counts )

matrixSize = data.frame( n_genes = nrow( counts ), n_cells = ncol( counts ), nnz = length( counts@x ) )
matrixSize
fwrite( matrixSize, file = paste0( tablePath, "scCDC_CorrectedCounts_Size.txt" ), sep = "\t" )

saveRDS( counts, file = paste0( dataPath, "NatComm_CorrectedCounts_scCDC.rds" ) )
