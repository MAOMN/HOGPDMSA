library(clusterProfiler)
library(org.Hs.eg.db)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_1/" )
tablePath = paste0( outputPath, "tables/" )
objectPath = paste0( outputPath, "objects/" )

writeTable <- function( x, file, quote = FALSE ) write.table( x, file = paste0( tablePath, file ), col.names = TRUE, row.names = FALSE, quote = quote, sep = "\t" )


# tested genes of 01_DifferentialExpression_Analysis.R (columns: gene, symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, significance), only protein-coding genes
deg = read.delim( paste0( tablePath, "DifferentialExpression/HOG_DEG_WithGeneSymbols.txt" ) )
deg <- deg[ deg$gene %in% readLines( paste0( dataPath, "ReferenceFiles/GENCODE_v47/GENCODE_v47_ProteinCodingGenes.txt" ) ), ]

# the gene set collections use Entrez gene ids
deg$entrez <- mapIds( org.Hs.eg.db, keys = deg$gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first" )[ deg$gene ]
deg <- deg[ !is.na( deg$entrez ), ]

# gene universe: all genes with an adjusted p-value, gene lists: differentially expressed genes
geneUniverse = unique( deg$entrez[ !is.na( deg$padj ) ] )
sig <- deg[ !is.na( deg$padj ) & deg$padj < 0.05, ]
geneSets = list( Up = unique( sig$entrez[ sig$log2FoldChange > 0 ] ),
                 Down = unique( sig$entrez[ sig$log2FoldChange < 0 ] ),
                 All = unique( sig$entrez ) )

# gene set collections: GO and Reactome from MSigDB 2026.1 (made by DataPreparation/ReferenceFiles/code/MSigDB_GeneSets.R) and KEGG
collections = c( "GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome" )
msigdb = sapply( setdiff( collections, "KEGG" ), function( db ) as.data.frame( readRDS( paste0( dataPath, "ReferenceFiles/GeneSets/MSigDB_2026.1_", db, ".rds" ) ) ), simplify = FALSE )
# KEGG pathways, made by DataPreparation/ReferenceFiles/code/KEGG_Pathways.R (KEGG has no version number, the article uses the download of 2026-09-05)
kegg = readRDS( paste0( dataPath, "ReferenceFiles/GeneSets/KEGG_hsa_Pathways.rds" ) )


# enrichment of a gene list in one gene set collection
runEnrichment <- function( genes, db )
{
    if( db == "KEGG" ) return( enrichKEGG( gene = genes, organism = kegg, keyType = "ncbi-geneid", universe = geneUniverse, minGSSize = 5 ) )

    m = msigdb[[ db ]]
    term2gene = unique( m[ , c( "gs_name", "ncbi_gene" ) ] )
    term2gene = term2gene[ !is.na( term2gene$ncbi_gene ), ]
    term2gene = data.frame( term = term2gene$gs_name, gene = as.character( term2gene$ncbi_gene ) )

    # GO terms only have the name with underscores, make them readable
    term2name = unique( m[ , c( "gs_name", "gs_description", "gs_subcollection" ) ] )
    goName = tolower( gsub( "_", " ", sub( "^(GOBP|GOCC|GOMF|KEGG|REACTOME|HP)_", "", term2name$gs_name ) ) )
    substr( goName, 1, 1 ) <- toupper( substr( goName, 1, 1 ) )
    term2name = data.frame( term = term2name$gs_name, name = ifelse( grepl( "^GO:", term2name$gs_subcollection ), goName, term2name$gs_description ) )

    enricher( gene = genes, universe = geneUniverse, minGSSize = 5, TERM2GENE = term2gene, TERM2NAME = term2name )
}

# table of enriched gene sets with gene symbols and link to the gene set
getResultTable <- function( res, db )
{
    if( is.null( res ) ) return( data.frame() )
    tab = as.data.frame( setReadable( res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID" ) )
    if( nrow( tab ) == 0 ) return( tab )

    if( db == "KEGG" )
    {
        tab$Description = sub( "\\s*-\\s*Homo sapiens\\s*\\(human\\)\\s*$", "", tab$Description )
        tab$gs_id = tab$ID
        tab$gs_exact_source = tab$ID
        tab$gs_url = paste0( "https://www.kegg.jp/pathway/", tab$ID )
    }
    else
    {
        m = unique( msigdb[[ db ]][ , c( "gs_name", "gs_id", "gs_exact_source", "gs_url" ) ] )
        cand = match( tab$ID, m$gs_name )
        tab$gs_id = m$gs_id[ cand ]
        tab$gs_exact_source = m$gs_exact_source[ cand ]
        tab$gs_url = m$gs_url[ cand ]
    }
    tab
}


# enrichment for up-regulated, down-regulated and all differentially expressed genes
overview = c()
resultTables = list()
for( direction in c( "Up", "Down", "All" ) )
{
    for( db in collections )
    {
        res = runEnrichment( genes = geneSets[[ direction ]], db = db )
        tab = getResultTable( res, db )
        writeTable( tab, paste0( "Enrichment/Enrichment_", db, "_", direction, "_DEG.txt" ) )
        if( direction == "All" ) resultTables[[ db ]] <- tab
        overview = rbind( overview, data.frame( db = db, direction = direction, n_tested = if( is.null( res ) ) 0 else nrow( res@result ),
                                                n_reported = nrow( tab ), n_sig = if( nrow( tab ) > 0 ) sum( tab$p.adjust < 0.05 ) else 0 ) )
    }
}
writeTable( overview, "Enrichment/Enrichment_Overview.txt" )

# enrichment of all differentially expressed genes, used by 04_Enrichment_Plots.R
saveRDS( resultTables, file = paste0( objectPath, "Enrichment_Results.rds" ) )
