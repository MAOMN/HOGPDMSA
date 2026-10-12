library(readxl)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# folder with the downloaded raw data (not part of this repository, see Data/RawDataInfo.txt)
rawDataPath <- "path/to/raw_data/"

# data folder
dataPath = paste0( home_directory, "Data/ExternalPDCohorts/" )


# donor table of the Lee data set (GEO series GSE148434), the input of the other scripts of the Lee data set
# gsm_id and sample_id: from the names of the twelve molecule files of GEO (<GSM>_<sample_id>_snRNA.molecule_info.h5)
# diagnosis, sex and age: from Table S1 of Lee et al., Sci Adv 9, eabo2467 (2023), not given in GEO
geoPath = paste0( rawDataPath, "ExternalPDCohorts/Lee_GSE148434/" )

files = sub( "_snRNA.molecule_info.h5$", "", list.files( geoPath, pattern = "_snRNA.molecule_info.h5$" ) )
stopifnot( length( files ) == 12 )
donors = data.table( gsm_id = sub( "_.*$", "", files ), sample_id = sub( "^[^_]*_", "", files ) )

# Table S1 (xlsx file of the supplementary material of the article, in the folder of the molecule files):
# one row for each donor, the header row has the column names
tableS1 = as.data.frame( read_excel( paste0( geoPath, "abo2467_Table_S1.xlsx" ), col_names = FALSE, col_types = "text", .name_repair = "minimal" ) )
header = which( tableS1[ , 2 ] == "Subject ID" )
stopifnot( length( header ) == 1 )
columns = match( c( "Subject ID", "Diagnosis", "Age of Death", "Sex" ), tableS1[ header, ] )
stopifnot( !anyNA( columns ) )
info = tableS1[ ( header + 1 ):nrow( tableS1 ), columns ]
colnames( info ) <- c( "sample_id", "diagnosis", "age", "sex" )
info = info[ info$sample_id %in% donors$sample_id, ]
stopifnot( nrow( info ) == 12, !anyDuplicated( info$sample_id ) )

donors = merge( donors, info, by = "sample_id" )
donors$diagnosis = c( Normal = "HC", Parkinson = "PD" )[ donors$diagnosis ]
donors$sex = c( Female = "F", Male = "M" )[ donors$sex ]
donors$age = as.integer( donors$age )
stopifnot( !anyNA( donors ) )

donors = donors[ order( gsm_id ), c( "gsm_id", "sample_id", "diagnosis", "sex", "age" ) ]
print( donors )

dir.create( dataPath, recursive = TRUE, showWarnings = FALSE )
fwrite( donors, file = paste0( dataPath, "Lee_Donors.txt" ), sep = "\t" )
