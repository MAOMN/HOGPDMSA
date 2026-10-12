library(fgsea)
library(BiocParallel)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_3/" )
tablePath = paste0( outputPath, "tables/" )
objectPath = paste0( outputPath, "objects/" )


states = c( "OPC", paste0( "OL", 1:6 ) )
comparisons = c( "MSA_vs_HC", "PD_vs_HC" )

# hallmark gene sets of MSigDB 2026.1, made by DataPreparation/ReferenceFiles/code/MSigDB_GeneSets.R
# columns: db, term_id, term_name, gene, source, source_collection, source_subcollection, source_version
hallmark = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneSets/MSigDB_2026.1_Hallmark.txt.gz" ) ) )


#
### GSEA for each cell state and comparison, genes are ranked by the Wald statistic of DESeq2
#

set.seed( 42 )
register( MulticoreParam( workers = 8 ) )

gseaTab = c()
universeTab = c()
for( comp in comparisons )
{
    for( state in states )
    {
        # results of 08_Pseudobulk_DESeq2_Analysis.R
        # columns: natural_label, contrast, case, reference, model_formula, gene, gene_type, is_protein_coding, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, ...
        deg = as.data.frame( fread( paste0( tablePath, "Pseudobulk_DESeq2/Pseudobulk_DESeq2_", comp, "_", state, ".txt.gz" ) ) )

        # protein-coding genes only
        rankTab = deg[ which( deg$is_protein_coding & !is.na( deg$padj ) & !is.na( deg$stat ) & deg$stat != 0 ), ]
        # sort by gene name first, so genes with exactly the same statistic always have the same order
        rankTab = rankTab[ order( rankTab$gene, method = "radix" ), ]
        stats = rankTab$stat
        names( stats ) = rankTab$gene
        stats = sort( stats, decreasing = TRUE )

        universeTab = rbind( universeTab, data.frame( natural_label = state, contrast = comp, n_input_rows = nrow( deg ),
                                                      n_protein_coding_padj_nonNA = sum( deg$is_protein_coding & !is.na( deg$padj ), na.rm = TRUE ),
                                                      n_analysis_universe = length( stats ), rank_metric = "DESeq2 Wald statistic" ) )

        # hallmark gene sets with 5 to 500 measured genes
        h = hallmark[ hallmark$gene %in% names( stats ), ]
        pathways = split( h$gene, h$term_id )
        pathways = pathways[ sapply( pathways, function( g ) length( unique( g ) ) ) >= 5 & sapply( pathways, function( g ) length( unique( g ) ) ) <= 500 ]

        res = as.data.frame( fgseaMultilevel( pathways = pathways, stats = stats, minSize = 5, maxSize = 500, eps = 0, nPermSimple = 10000 ) )
        res$leadingEdge = sapply( res$leadingEdge, paste, collapse = ";" )
        colnames( res )[ colnames( res ) == "pathway" ] = "term_id"
        res$term_name = hallmark$term_name[ match( res$term_id, hallmark$term_id ) ]
        gseaTab = rbind( gseaTab, data.frame( natural_label = state, contrast = comp, db = "Hallmark", res[ , c( "term_id", "term_name", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge" ) ],
                                              rank_metric = "DESeq2 Wald statistic", minGSSize = 5, maxGSSize = 500, nPermSimple = 10000 ) )
    }
}

gseaTab$significant <- !is.na( gseaTab$padj ) & gseaTab$padj < 0.05
gseaTab$neg_log10_fdr <- -log10( pmax( gseaTab$padj, .Machine$double.xmin ) )

# readable names: HALLMARK_TNFA_SIGNALING_VIA_NFKB -> Tnfa signaling via nfkb
hallmarkName <- tolower( gsub( "_", " ", sub( "^HALLMARK_", "", gseaTab$term_name ) ) )
substr( hallmarkName, 1, 1 ) <- toupper( substr( hallmarkName, 1, 1 ) )
gseaTab$display_term <- hallmarkName

fwrite( gseaTab, file = paste0( tablePath, "GSEA_Hallmark/GSEA_Hallmark_AllResults.txt.gz" ), sep = "\t" )
fwrite( universeTab, file = paste0( tablePath, "GSEA_Hallmark/GSEA_Hallmark_RankedGenes_Overview.txt" ), sep = "\t" )

# number of significant gene sets
overview = c()
for( comp in comparisons )
{
    for( state in states )
    {
        cand = which( gseaTab$natural_label == state & gseaTab$contrast == comp )
        overview = rbind( overview, data.frame( natural_label = state, contrast = comp, n_results = length( cand ), n_fdr_significant = sum( gseaTab$significant[ cand ] ) ) )
    }
}
fwrite( overview, file = paste0( tablePath, "GSEA_Hallmark/GSEA_Hallmark_Overview.txt" ), sep = "\t" )

# the 4 most significant gene sets of each cell state and comparison
sig = gseaTab[ gseaTab$significant, ]
sig = sig[ order( sig$natural_label, sig$contrast, sig$padj, -abs( sig$NES ), -abs( sig$ES ), -sig$size, method = "radix" ), ]
groupId = paste( sig$natural_label, sig$contrast )
topTerms = sig[ unlist( lapply( unique( groupId ), function( g ) which( groupId == g )[ 1:min( 4, sum( groupId == g ) ) ] ) ), ]
fwrite( topTerms, file = paste0( tablePath, "GSEA_Hallmark/GSEA_Hallmark_Top4_PerStateAndComparison.txt" ), sep = "\t" )

# data of the plot in 13_GSEA_Hallmark_Plots.R
saveRDS( list( gseaTab = gseaTab, topTerms = topTerms ), file = paste0( objectPath, "GSEA_Hallmark_PlotData.rds" ) )
