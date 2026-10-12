library(clusterProfiler)
library(gson)

# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/ReferenceFiles/GeneSets/" )


# KEGG has no version number: the article uses the download of 2026-09-05, a new download can give slightly different results for Figure 1
# the KEGG API is for academic use only (https://www.kegg.jp/kegg/rest/)
# list with: KEGGPATHID2EXTID (pathway id - NCBI gene id), KEGGPATHID2NAME (pathway id - pathway name)
keggDownload <- download_KEGG( "hsa" )

# gene set object that enrichKEGG accepts instead of downloading the pathways again
pathwayGenes = keggDownload$KEGGPATHID2EXTID
colnames( pathwayGenes ) <- c( "gsid", "gene" )
pathwayNames = keggDownload$KEGGPATHID2NAME
colnames( pathwayNames ) <- c( "gsid", "name" )

kegg <- gson( gsid2gene = pathwayGenes, gsid2name = pathwayNames, species = "hsa", gsname = "KEGG", keytype = "ncbi-geneid", accessed_date = as.character( Sys.Date() ) )

dir.create( dataPath, recursive = TRUE, showWarnings = FALSE )
saveRDS( kegg, file = paste0( dataPath, "KEGG_hsa_Pathways.rds" ) )
