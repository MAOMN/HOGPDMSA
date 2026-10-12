library(dynamicTreeCut)
library(fgsea)
library(clusterProfiler)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/" )

# results of the HOG cells: ranking of the nuclear genes by the Wald statistic (Figure 2)
hogPath = paste0( home_directory, "Figure_2/" )

# output folder
outputPath = paste0( home_directory, "Figure_5/" )
tablePath = paste0( outputPath, "tables/Modules/" )
objectPath = paste0( outputPath, "objects/" )


units = c( "MSA_L1", "MSA_L3", "PD_L1", "PD_L3" )

# fitted curves of the genes that differ between disease and HC (made by 02_Disease_Dynamics_ConditionTest_Analysis.R)
# list with one entry per disease and lineage: genes, raw_response, log1p_response, scaled (z-score of log1p, genes x 100: 50 cell groups HC, 50 cell groups disease),
# time (median pseudotime of the cell groups, 0 = start and 1 = end of the compared range), median_pseudotime
curves <- readRDS( paste0( objectPath, "Disease_Dynamics_Curves.rds" ) )

# results of the condition test (made by 02_Disease_Dynamics_ConditionTest_Analysis.R)
condTab <- readRDS( paste0( objectPath, "Disease_Dynamics_ConditionTest.rds" ) )

clean_symbol <- function( x )
{
    x = toupper( trimws( as.character( x ) ) )
    x[ is.na( x ) ] = ""
    x
}


#
### Gene modules: hierarchical clustering of the z-score curves (HC and disease curve together), modules by dynamic tree cut
#

assignTab = c()
treeTab = c()
mainCurveTab = c()
stageTab = c()
for( u in units )
{
    x = curves[[ u ]]
    d = dist( x$scaled )
    tree = hclust( d, method = "ward.D2" )
    labels = cutreeDynamic( dendro = tree, distM = as.matrix( d ), deepSplit = 0, minClusterSize = 20 )

    assigned = labels > 0
    module = ifelse( assigned, sprintf( "D%02d", labels ), "Unassigned" )
    a = data.frame( unit_id = u, disease = x$disease, reference = "HC", lineage_id = x$lineage_id, endpoint = x$endpoint, gene = x$genes,
                    raw_cluster = as.integer( labels ), assigned = assigned, module = module, module_id = paste( u, module, sep = "__" ) )
    assignTab = rbind( assignTab, a )
    treeTab = rbind( treeTab, data.frame( unit_id = u, tree_order = 1:length( x$genes ), gene = x$genes[ tree$order ] ) )

    #mean curve of each module, separately for HC (points 1-50) and disease (points 51-100)
    for( mid in unique( a$module_id[ a$assigned ] ) )
    {
        members = a$gene[ a$module_id == mid ]
        for( ci in 1:2 )
        {
            cols = ( ci - 1 ) * 50 + 1:50
            mainCurveTab = rbind( mainCurveTab, data.frame( unit_id = u, module_id = mid, condition = c( "HC", x$disease )[ ci ], normalized_pseudotime = x$time, genes = length( members ),
                                                            mean_scaled = colMeans( x$scaled[ members, cols, drop = FALSE ] ),
                                                            median_scaled = apply( x$scaled[ members, cols, drop = FALSE ], 2, median ),
                                                            mean_response = colMeans( x$raw_response[ members, cols, drop = FALSE ] ),
                                                            median_response = apply( x$raw_response[ members, cols, drop = FALSE ], 2, median ) ) )
        }
    }

    sizes = table( labels[ assigned ] )
    stageTab = rbind( stageTab, data.frame( unit_id = u, input_genes = length( x$genes ), assigned_genes = sum( assigned ), unassigned_genes = sum( !assigned ), modules = length( sizes ),
                                            min_module_size = min( sizes ), median_module_size = median( as.numeric( sizes ) ), max_module_size = max( sizes ) ) )
}

fwrite( assignTab, file = paste0( tablePath, "Disease_Dynamics_Modules_Genes.txt.gz" ), sep = "\t" )
fwrite( treeTab, file = paste0( tablePath, "Disease_Dynamics_Modules_TreeOrder.txt.gz" ), sep = "\t" )
fwrite( mainCurveTab, file = paste0( tablePath, "Disease_Dynamics_Modules_MeanCurves.txt.gz" ), sep = "\t" )
fwrite( stageTab, file = paste0( tablePath, "Disease_Dynamics_Modules_Sizes.txt.gz" ), sep = "\t" )

assignedTab = assignTab[ assignTab$assigned, ]


#
### Are the genes of a module changed in the same direction in the HOG cells treated with alpha-syn PFF? (GSEA of the modules in the HOG gene ranking)
#

# ranking of the nuclear encoded protein-coding genes of the HOG cells by Wald statistic (made by 03_GSEA_MitoCarta_Analysis.R of Figure 2)
# columns: rank, gene, stat, stat_plot
hogRank = as.data.frame( fread( paste0( hogPath, "tables/GSEA_OXPHOS/GSEA_NuclearGenes_WaldRanking.txt.gz" ) ) )
hogRank$symbol = clean_symbol( hogRank$gene )
hogRank = hogRank[ hogRank$symbol != "" & is.finite( hogRank$stat ) & !duplicated( hogRank$symbol ), ]
hogStats = hogRank$stat
names( hogStats ) = hogRank$symbol
hogStats = sort( hogStats, decreasing = TRUE )

hogTab = c()
testableTab = c()
for( u in units )
{
    a = assignedTab[ assignedTab$unit_id == u, ]
    sets = lapply( split( clean_symbol( a$gene ), a$module_id ), unique )

    #a module needs between 15 and 500 genes in the HOG ranking
    overlap = sapply( sets, function( g ) length( intersect( g, names( hogStats ) ) ) )
    testableTab = rbind( testableTab, data.frame( unit_id = u, module_id = names( sets ), module_genes = sapply( sets, length ), HOG_overlap = overlap, testable = overlap >= 15 & overlap <= 500 ) )

    set.seed( 42 )
    res = fgseaMultilevel( pathways = sets, stats = hogStats, minSize = 15, maxSize = 500, eps = 0, nproc = 1, nPermSimple = 10000 )
    res = as.data.frame( res )
    res$leadingEdge = sapply( res$leadingEdge, paste, collapse = ";" )
    colnames( res )[ colnames( res ) == "pathway" ] = "module_id"
    colnames( res )[ colnames( res ) == "padj" ] = "padj_fgsea_call"
    # p-values adjusted over the modules of one disease and lineage
    res$padj_unit_BH = p.adjust( res$pval, method = "BH" )
    res$unit_id = u
    res$correction_family = paste0( u, "__all_testable_DTC_modules" )
    res$HOG_FDR05 = is.finite( res$padj_unit_BH ) & res$padj_unit_BH < 0.05
    hogTab = rbind( hogTab, res )
}

fwrite( hogTab, file = paste0( tablePath, "Disease_Dynamics_Modules_HOG_GSEA.txt.gz" ), sep = "\t", na = "NA" )
fwrite( testableTab, file = paste0( tablePath, "Disease_Dynamics_Modules_HOG_GenesInRanking.txt.gz" ), sep = "\t", na = "NA" )


#
### GO biological process enrichment of each module
#

# GO BP gene sets of MSigDB 2026.1, made by DataPreparation/ReferenceFiles/code/MSigDB_GeneSets.R
gobp = as.data.frame( readRDS( paste0( dataPath, "ReferenceFiles/GeneSets/MSigDB_2026.1_GO_BP.rds" ) ) )
term2gene = unique( data.frame( term_id = gobp$gs_id, gene = clean_symbol( gobp$gene_symbol ) ) )
term2name = unique( data.frame( term_id = gobp$gs_id, term_name = gobp$gs_name ) )

# gene universe: all genes with a valid test for this disease and lineage (that are in at least one GO BP gene set)
goTab = c()
for( mid in unique( assignedTab$module_id ) )
{
    a = assignedTab[ assignedTab$module_id == mid, ]
    u = a$unit_id[ 1 ]
    universeAll = unique( clean_symbol( condTab$gene[ condTab$unit_id == u & condTab$finite_test ] ) )
    universe = intersect( universeAll, unique( term2gene$gene ) )
    genesAll = unique( clean_symbol( a$gene ) )
    genes = intersect( genesAll, universe )

    res = enricher( gene = genes, universe = universe, TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 15, maxGSSize = 500 )
    res = as.data.frame( res )

    colnames( res )[ match( c( "ID", "Description", "p.adjust", "Count" ), colnames( res ) ) ] = c( "term_id", "term_name", "padj_clusterProfiler", "overlap_n" )
    res$padj_within_module_BH = p.adjust( res$pvalue, method = "BH" )
    res$module_id = mid
    res$unit_id = u
    res$disease = a$disease[ 1 ]
    res$lineage_id = a$lineage_id[ 1 ]
    res$module_genes = length( genesAll )
    res$query_genes_mapped = length( genes )
    res$background_genes = length( universeAll )
    res$background_mapped = length( universe )
    res$correction_family = paste0( mid, "__all_overlapping_GO_BP_terms" )
    res$overlap_genes = gsub( "/", ";", res$geneID, fixed = TRUE )
    res$FDR05 = is.finite( res$padj_within_module_BH ) & res$padj_within_module_BH < 0.05
    goTab = rbind( goTab, res )
}
rownames( goTab ) <- NULL

fwrite( goTab, file = paste0( tablePath, "Disease_Dynamics_Modules_GO_BP.txt.gz" ), sep = "\t", na = "NA" )

# overview of the modules
moduleTab = unique( assignedTab[ , c( "unit_id", "module_id", "module" ) ] )
moduleTab$genes = sapply( moduleTab$module_id, function( m ) sum( assignedTab$module_id == m ) )
cand = match( moduleTab$module_id, testableTab$module_id )
moduleTab$HOG_overlap = testableTab$HOG_overlap[ cand ]
moduleTab$HOG_testable = testableTab$testable[ cand ]
cand = match( moduleTab$module_id, hogTab$module_id )
moduleTab$HOG_NES = hogTab$NES[ cand ]
moduleTab$HOG_pvalue = hogTab$pval[ cand ]
moduleTab$HOG_FDR = hogTab$padj_unit_BH[ cand ]
moduleTab$HOG_FDR05 = hogTab$HOG_FDR05[ cand ]
moduleTab$GO_BP_tested_terms = sapply( moduleTab$module_id, function( m ) sum( goTab$module_id == m ) )
moduleTab$GO_BP_FDR05_terms = sapply( moduleTab$module_id, function( m ) sum( goTab$FDR05[ goTab$module_id == m ] ) )
moduleTab = moduleTab[ order( match( moduleTab$unit_id, units ), moduleTab$module ), ]

# names of the modules in the figures: C1, C2, ... within each disease and lineage
moduleTab$unit_order = match( moduleTab$unit_id, units )
moduleTab$module_order = as.integer( sub( "^D", "", moduleTab$module ) )
moduleTab$display_cluster = unlist( lapply( units, function( u ) paste0( "C", 1:sum( moduleTab$unit_id == u ) ) ) )
moduleTab$panel_order = 1:nrow( moduleTab )
moduleTab$display_unit = gsub( "_", " ", moduleTab$unit_id )

fwrite( moduleTab, file = paste0( tablePath, "Disease_Dynamics_Modules_Overview.txt.gz" ), sep = "\t", na = "NA" )

# the 6 most significant GO terms of each module
topGoTab = c()
for( mid in sort( unique( goTab$module_id[ goTab$FDR05 ] ) ) )
{
    g = goTab[ goTab$module_id == mid & goTab$FDR05, ]
    g = g[ order( g$padj_within_module_BH, g$pvalue ), ]
    g = g[ 1:min( 6, nrow( g ) ), ]
    g$display_cluster = moduleTab$display_cluster[ moduleTab$module_id == mid ]
    topGoTab = rbind( topGoTab, g )
}
fwrite( topGoTab, file = paste0( tablePath, "Disease_Dynamics_Modules_GO_BP_Top6.txt" ), sep = "\t", na = "NA" )


#
### Results for the plots (04_Disease_Dynamics_Modules_Plots.R)
#

saveRDS( list( assignTab = assignTab, moduleTab = moduleTab, goTab = goTab ), file = paste0( objectPath, "Disease_Dynamics_Modules_PlotData.rds" ) )
