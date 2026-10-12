library(data.table)
library(openxlsx)

# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/" )

# results of the comparison of disease and HC along the trajectories (Figure 5)
dynPath = paste0( home_directory, "Figure_5/" )

# output folder
outputPath = paste0( home_directory, "SupplementaryTables/Table_S4/" )
tablePath = paste0( outputPath, "tables/" )


#
### Supplementary Table 4 of the manuscript: FDR of selected GO biological processes in the modules, one sheet (only a new layout of the results in tables/)
#
# rows: 5 respiration and 5 protein folding processes, then the 6 processes with the lowest FDR of every module that is not C1 (processes already listed are not repeated)
# columns: the modules with at least one process with FDR < 0.05
# cells: the FDR of the process in the module, also when it is not significant

# GO BP enrichment of all modules (made by 01_GO_BP_Enrichment_DiseaseModules_Table.R)
go = fread( paste0( tablePath, "Supplementary_Table_S4_GO_BP_enrichment.txt.gz" ) )
go[ , module := paste( disease, lineage, module_label ) ]
go[ , unit_id := paste( disease, lineage, sep = "_" ) ]
stopifnot( !anyDuplicated( go[ , .( module, msigdb_gene_set_name ) ] ), all( is.finite( go$fdr_bh_within_module ) ) )

modOrder = unique( go[ , .( disease, lineage, module_label, module ) ] )
modOrder = modOrder[ order( match( disease, c( "MSA", "PD" ) ), lineage, module_label ) ]
sigGo = go[ fdr_bh_within_module < 0.05 ]
modules = modOrder$module[ modOrder$module %in% sigGo$module ]

respiration = c( "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_ELECTRON_TRANSPORT_CHAIN", "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT", "GOBP_CELLULAR_RESPIRATION",
                 "GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS" )
folding = c( "GOBP_PROTEIN_FOLDING", "GOBP_PROTEIN_REFOLDING", "GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN", "GOBP_RESPONSE_TO_HEAT", "GOBP_CELLULAR_RESPONSE_TO_HEAT" )
chosen = c( respiration, folding )

otherTop = c()
for( m in modules[ !grepl( " C1$", modules ) ] )
{
    x = sigGo[ module == m ][ order( fdr_bh_within_module, p_value, msigdb_gene_set_name ) ]
    otherTop = c( otherTop, head( x$msigdb_gene_set_name, 6 ) )
}
otherTop = otherTop[ !duplicated( otherTop ) & !otherTop %in% chosen ]
processes = c( chosen, otherTop )
stopifnot( all( processes %in% go$msigdb_gene_set_name ) )

# the enrichment returns only processes with a gene in the module; a process with 15 to 500 genes in the gene universe and none in the module has FDR = 1
# (gene universe as in 03_Disease_Dynamics_Modules_Analysis.R of Figure 5)
clean_symbol <- function( x )
{
    x = toupper( trimws( as.character( x ) ) )
    x[ is.na( x ) ] = ""
    x
}
# condition tests (made by 02_Disease_Dynamics_ConditionTest_Analysis.R of Figure 5)
condTab = readRDS( paste0( dynPath, "objects/Disease_Dynamics_ConditionTest.rds" ) )
# GO BP gene sets of MSigDB 2026.1, made by DataPreparation/ReferenceFiles/code/MSigDB_GeneSets.R
gobp = as.data.frame( readRDS( paste0( dataPath, "ReferenceFiles/GeneSets/MSigDB_2026.1_GO_BP.rds" ) ) )
term2gene = unique( data.frame( term_name = gobp$gs_name, gene = clean_symbol( gobp$gene_symbol ) ) )

sizes = c()
for( u in unique( go$unit_id ) )
{
    universe = intersect( unique( clean_symbol( condTab$gene[ condTab$unit_id == u & condTab$finite_test ] ) ), unique( term2gene$gene ) )
    t2g = term2gene[ term2gene$gene %in% universe, ]
    sizes = rbind( sizes, data.table( unit_id = u, msigdb_gene_set_name = names( table( t2g$term_name ) ), size = as.integer( table( t2g$term_name ) ) ) )
}
# the sizes must be the sizes of the enrichment results
check = merge( go[ , .( unit_id, msigdb_gene_set_name, gene_set_size_in_universe ) ], sizes, by = c( "unit_id", "msigdb_gene_set_name" ) )
stopifnot( nrow( check ) == nrow( go ), all( check$size == check$gene_set_size_in_universe ) )

# readable names
label = function( x )
{
    x = tolower( sub( "^GOBP_", "", x ) )
    x = gsub( "_", " ", x )
    x = gsub( "\\bt cell", "T cell", x )
    for( w in c( "atp", "nadh", "dna", "rna", "ii" ) ) x = gsub( paste0( "\\b", w, "\\b" ), if( w == "ii" ) "II" else toupper( w ), x )
    paste0( toupper( substr( x, 1, 1 ) ), substr( x, 2, nchar( x ) ) )
}

tab = data.frame( Biological_Process = label( processes ) )
for( m in modules )
{
    fdr = go[ module == m ][ match( processes, msigdb_gene_set_name ) ]$fdr_bh_within_module
    s = sizes[ unit_id == go[ module == m ]$unit_id[ 1 ] ][ match( processes, msigdb_gene_set_name ) ]$size
    noGene = is.na( fdr ) & !is.na( s ) & s >= 15 & s <= 500
    fdr[ noGene ] = 1
    # processes with fewer than 15 or more than 500 genes in the gene universe are not tested: the cell says Not tested
    stopifnot( identical( is.na( fdr ), is.na( s ) | s < 15 | s > 500 ) )
    tab[[ m ]] = fdr
}

# layout as Supplementary Table 1: header in the first row in bold, thin border around every cell, Arial 10
wb <- createWorkbook( creator = "" )
modifyBaseFont( wb, fontSize = 10, fontName = "Arial" )
sheetName = "Supplementary Table 4"
addWorksheet( wb, sheetName )
# column names with the lineage written out (MSA L1 C2 -> MSA Lineage 1 C2)
names( tab )[ -1 ] = sub( " L([0-9]) ", " Lineage \\1 ", names( tab )[ -1 ] )
writeData( wb, sheetName, tab, startRow = 1, startCol = 1 )
# processes that are not tested in a module
for( k in 2:ncol( tab ) ) for( r in which( is.na( tab[[ k ]] ) ) ) writeData( wb, sheetName, "Not tested", startCol = k, startRow = r + 1 )

cellStyle = function( ... ) createStyle( fontName = "Arial", fontSize = 10, border = "TopBottomLeftRight", borderColour = "black", borderStyle = "thin", ... )
rows = 2:( nrow( tab ) + 1 )
addStyle( wb, sheetName, cellStyle( textDecoration = "bold" ), rows = 1, cols = 1:ncol( tab ), gridExpand = TRUE )
addStyle( wb, sheetName, cellStyle(), rows = rows, cols = 1, gridExpand = TRUE )
addStyle( wb, sheetName, cellStyle( numFmt = "0.00E+00" ), rows = rows, cols = 2:ncol( tab ), gridExpand = TRUE )
# openxlsx stores each width 0.71 larger (52 becomes 52.71)
setColWidths( wb, sheetName, cols = 1:ncol( tab ), widths = c( 52, rep( 17, length( modules ) ) ) )
# row height of Arial 10 (openxlsx writes the default height of Calibri 11)
wb$worksheets[[ 1 ]]$sheetFormatPr <- '<sheetFormatPr defaultColWidth="11.42578125" defaultRowHeight="12.75"/>'
# remove links to drawing files that openxlsx does not write (macOS Quick Look cannot open the file otherwise)
wb$worksheets_rels[[ 1 ]] <- wb$worksheets_rels[[ 1 ]][ !grepl( "drawing", wb$worksheets_rels[[ 1 ]], ignore.case = TRUE ) ]

saveWorkbook( wb, paste0( outputPath, "Supplementary_Table_4.xlsx" ), overwrite = TRUE )
