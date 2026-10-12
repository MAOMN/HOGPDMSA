library(data.table)
library(openxlsx)

# change to your project directory
home_directory <- "path/to/project_directory/"

# output folder
outputPath = paste0( home_directory, "SupplementaryTables/Table_S1/" )
tablePath = paste0( outputPath, "tables/" )


#
### Supplementary Table 1 of the manuscript: the results of the selected TFs with FDR < 0.05, one sheet (only a new layout of the results in tables/)
#

# all TF activities (made by 01_TF_Activity_ExternalPDCohorts_Analysis.R) and the selected TFs (made by 02_TF_Selection_Analysis.R)
tf = fread( paste0( tablePath, "Supplementary_Table_S1_TF_activity.txt" ) )
selected = fread( paste0( tablePath, "ExternalPDCohorts_SelectedTFs.txt" ) )
datasets = c( "Martirosyan", "Smajic", "Lee" )
pops = c( "OPC", "OL1", "OL2", "OL3", "OL4", "OL5" )

# one row for each selected TF and population with FDR < 0.05, ordered by data set, order of the selection and population
d = merge( selected[ , .( dataset, TF, display_order ) ], tf, by = c( "dataset", "TF" ) )
stopifnot( all( d$test_status == "TESTED" ) )
d = d[ is.finite( FDR_BH_within_population ) & FDR_BH_within_population < 0.05 ]
d = d[ order( match( dataset, datasets ), display_order, match( population, pops ) ) ]

tab = data.frame( Dataset = ifelse( d$dataset == "Smajic", "Smajić", d$dataset ), TF = d$TF, Population = d$population, Activity = d$VIPER_NES, FDR = d$FDR_BH_within_population )

# layout: header in the first row in bold, thin border around every cell, Arial 10
wb <- createWorkbook( creator = "" )
modifyBaseFont( wb, fontSize = 10, fontName = "Arial" )
sheetName = "Supplementary Table 1"
addWorksheet( wb, sheetName )
writeData( wb, sheetName, tab, startRow = 1, startCol = 1 )

cellStyle = function( ... ) createStyle( fontName = "Arial", fontSize = 10, border = "TopBottomLeftRight", borderColour = "black", borderStyle = "thin", ... )
rows = 2:( nrow( tab ) + 1 )
addStyle( wb, sheetName, cellStyle( textDecoration = "bold" ), rows = 1, cols = 1:ncol( tab ), gridExpand = TRUE )
addStyle( wb, sheetName, cellStyle(), rows = rows, cols = 1:3, gridExpand = TRUE )
addStyle( wb, sheetName, cellStyle( numFmt = "0.00" ), rows = rows, cols = 4, gridExpand = TRUE )
addStyle( wb, sheetName, cellStyle( numFmt = "0.00E+00" ), rows = rows, cols = 5, gridExpand = TRUE )
# openxlsx stores each width 0.71 larger (14 becomes 14.71)
setColWidths( wb, sheetName, cols = 1:ncol( tab ), widths = c( 14, 9, 12, 9, 11 ) )
# row height of Arial 10 (openxlsx writes the default height of Calibri 11)
wb$worksheets[[ 1 ]]$sheetFormatPr <- '<sheetFormatPr defaultColWidth="11.42578125" defaultRowHeight="12.75"/>'
# remove links to drawing files that openxlsx does not write (macOS Quick Look cannot open the file otherwise)
wb$worksheets_rels[[ 1 ]] <- wb$worksheets_rels[[ 1 ]][ !grepl( "drawing", wb$worksheets_rels[[ 1 ]], ignore.case = TRUE ) ]

saveWorkbook( wb, paste0( outputPath, "Supplementary_Table_1.xlsx" ), overwrite = TRUE )
