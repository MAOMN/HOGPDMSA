library(msigdbr)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/ReferenceFiles/GeneSets/" )


# gene set collections of MSigDB 2026.1 (human) from the package msigdbr 26.1.0, which downloads them once from https://zenodo.org/records/18968178
dir.create( dataPath, recursive = TRUE, showWarnings = FALSE )

# GO (biological process, cellular component, molecular function) and Reactome are saved as returned by msigdbr
collections = data.frame( name = c( "GO_BP", "GO_CC", "GO_MF", "Reactome" ), collection = c( "C5", "C5", "C5", "C2" ),
                          subcollection = c( "GO:BP", "GO:CC", "GO:MF", "CP:REACTOME" ) )
for( i in 1:nrow( collections ) )
{
    geneSets <- msigdbr( species = "Homo sapiens", collection = collections$collection[ i ], subcollection = collections$subcollection[ i ] )
    print( c( collections$name[ i ], unique( geneSets$db_version ), length( unique( geneSets$gs_name ) ) ) )
    saveRDS( geneSets, file = paste0( dataPath, "MSigDB_2026.1_", collections$name[ i ], ".rds" ) )
}

# hallmark gene sets as table with one row per gene set and gene (gene symbols)
hallmark <- as.data.table( msigdbr( species = "Homo sapiens", collection = "H" ) )
hallmarkTab <- unique( hallmark[ , .( db = "Hallmark", term_id = gs_id, term_name = gs_name, gene = gene_symbol, source = "msigdbr",
                                      source_collection = "H", source_subcollection = "", source_version = db_version ) ] )
hallmarkTab <- hallmarkTab[ !is.na( gene ) & gene != "" & !is.na( term_id ) & term_id != "" ]

fwrite( hallmarkTab, file = paste0( dataPath, "MSigDB_2026.1_Hallmark.txt.gz" ), sep = "\t" )
