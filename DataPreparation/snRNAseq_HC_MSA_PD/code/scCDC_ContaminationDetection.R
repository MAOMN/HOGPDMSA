library(scCDC)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# folder with the downloaded raw data (not part of this repository, see Data/snRNAseq_HC_MSA_PD/README.txt)
rawDataPath <- "path/to/raw_data/"

# Python with the package anndata, needed to read the h5ad file (change to your Python, see PythonPackages.txt)
pythonPath <- "path/to/python"

# output folder
outputPath = paste0( home_directory, "DataPreparation/snRNAseq_HC_MSA_PD/" )
codePath = paste0( outputPath, "code/" )
tablePath = paste0( outputPath, "tables/" )

set.seed( 42 )


# read the snRNA-seq data set of the authors
# gives: obs (cell annotation with donor and diagnosis), genes, function getDonorObject
source( paste0( codePath, "readAuthorCounts.R" ) )
length( genes )
table( obs$diagnosis )


#
### Detection of contaminating genes (ambient RNA) with scCDC, separately for each donor
#

donorTab = unique( obs[ , c( "diagnosis", "donor" ) ] )
donorTab = donorTab[ order( donorTab$diagnosis, donorTab$donor, method = "radix" ), ]
nrow( donorTab )

geneTab = list()
overview = c()
for( k in 1:nrow( donorTab ) )
{
    donor = donorTab$donor[ k ]
    diagnosis = donorTab$diagnosis[ k ]
    obj <- getDonorObject( donor )
    print( paste( k, donor, ncol( obj ), "cells" ) )

    #contaminating genes (GCGs): genes that are found in all cell types although they should be restricted to some
    detection <- ContaminationDetection( seuratobject = obj, restriction_factor = 0.5, sample_name = paste( "NatComm", diagnosis, donor, sep = "__" ), min.cell = 100, percent.cutoff = 0.2,
                                         out_path.plot = NULL, out_path.table = NULL, only.cont_genes = TRUE )

    #contamination ratio of the donor
    ratio <- ContaminationQuantification( object = obj, cont_genes = rownames( detection ), auc_thres = 0.9, min.cell = 100, top10_gcg = FALSE, layer = "counts" )

    #one column per cell type of the donor and the mean distance
    geneTab[[ k ]] = data.frame( diagnosis = diagnosis, sample_id = donor, gene = rownames( detection ), detection,
                                 n_cells = ncol( obj ), n_features = nrow( obj ), n_clusters = length( levels( Idents( obj ) ) ), check.names = FALSE, row.names = NULL )
    overview = rbind( overview, data.frame( diagnosis = diagnosis, sample_id = donor, n_cells = ncol( obj ), n_features = nrow( obj ), n_clusters = length( levels( Idents( obj ) ) ),
                                            n_gcg = nrow( detection ), contamination_ratio = as.numeric( ratio ), gcg_genes = paste( rownames( detection ), collapse = ";" ) ) )

    #save the tables after each donor
    fwrite( rbindlist( geneTab, fill = TRUE ), file = paste0( tablePath, "scCDC_ContaminatingGenes.txt" ), sep = "\t" )
    fwrite( overview, file = paste0( tablePath, "scCDC_Detection_Overview.txt" ), sep = "\t" )

    rm( obj )
    gc()
}
overview[ , 1:7 ]
