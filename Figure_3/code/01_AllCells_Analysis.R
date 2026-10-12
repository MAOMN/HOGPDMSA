library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_3/" )
tablePath = paste0( outputPath, "tables/AllCells/" )


# cell annotation of the snRNA-seq data set (HC, MSA and PD donors), cell type labels and UMAP coordinates are the ones of the authors of the data set
# made by DataPreparation/snRNAseq_HC_MSA_PD/code/CellAnnotation_AuthorLabels_UMAP.py from the data of Rydbirk et al. (2026)
# columns: cell_id, sample_id, diagnosis, labels_major, labels_minor, UMAP_1, UMAP_2
cells = as.data.frame( fread( paste0( dataPath, "snRNAseq_HC_MSA_PD/NatComm_CellAnnotation.txt.gz" ) ) )

# number of cells per cell type and per diagnosis
for( field in c( "labels_major", "labels_minor", "diagnosis" ) )
{
    counts = as.data.frame( table( cells[ , field ] ) )
    colnames( counts ) <- c( field, "n_cells" )
    counts = counts[ order( -counts$n_cells ), ]
    write.table( counts, file = paste0( tablePath, "AllCells_CellNumbers_", field, ".txt" ), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t" )
}
