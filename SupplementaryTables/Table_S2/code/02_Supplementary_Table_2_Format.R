library(data.table)
library(openxlsx)

# change to your project directory
home_directory <- "path/to/project_directory/"


# output folder
outputPath = paste0( home_directory, "SupplementaryTables/Table_S2/" )
tablePath = paste0( outputPath, "tables/" )


#
### Supplementary Table 2 of the manuscript: all genes of tables/, one sheet (only a new layout of the results in tables/)
#

writeTable <- function( tab, number, numericCols, widths )
{
    wb <- createWorkbook( creator = "" )
    modifyBaseFont( wb, fontSize = 10, fontName = "Arial" )
    sheetName = paste( "Supplementary Table", number )
    addWorksheet( wb, sheetName )
    writeData( wb, sheetName, tab, startRow = 1, startCol = 1, keepNA = FALSE )
    cellStyle = function( ... ) createStyle( fontName = "Arial", fontSize = 10, border = "TopBottomLeftRight", borderColour = "black", borderStyle = "thin", ... )
    rows = 2:( nrow( tab ) + 1 )
    addStyle( wb, sheetName, cellStyle( textDecoration = "bold" ), rows = 1, cols = 1:ncol( tab ), gridExpand = TRUE )
    addStyle( wb, sheetName, cellStyle(), rows = rows, cols = setdiff( 1:ncol( tab ), numericCols ), gridExpand = TRUE )
    addStyle( wb, sheetName, cellStyle( numFmt = "0.00E+00" ), rows = rows, cols = numericCols, gridExpand = TRUE )
    setColWidths( wb, sheetName, cols = 1:ncol( tab ), widths = widths )
    # row height of Arial 10 (openxlsx writes the default height of Calibri 11)
    wb$worksheets[[ 1 ]]$sheetFormatPr <- '<sheetFormatPr defaultColWidth="11.42578125" defaultRowHeight="12.75"/>'
    # openxlsx links the sheet to drawing files that it does not write; Excel ignores the links, but the preview of macOS (Quick Look) fails
    wb$worksheets_rels[[ 1 ]] <- wb$worksheets_rels[[ 1 ]][ !grepl( "drawing", wb$worksheets_rels[[ 1 ]], ignore.case = TRUE ) ]
    saveWorkbook( wb, paste0( outputPath, "Supplementary_Table_", number, ".xlsx" ), overwrite = TRUE )
}

x = as.data.frame( fread( paste0( tablePath, "Supplementary_Table_S2_HC_dynamic_genes.txt.gz" ) ) )
tab = data.frame( Lineage = sub( "^L", "Lineage ", x$lineage ), Gene = x$gene, FDR = x$fdr, Module = x$module, check.names = FALSE )
writeTable( tab, 2, numericCols = 3, widths = c( 10, 16, 10, 8 ) )
nrow( tab )
