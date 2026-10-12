library(fgsea)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
degPath = paste0( home_directory, "Figure_1/" )
outputPath = paste0( home_directory, "Figure_2/" )
tablePath = paste0( outputPath, "tables/GSEA_OXPHOS/" )
objectPath = paste0( outputPath, "objects/" )

writeTable <- function( x, file ) fwrite( x, file = paste0( tablePath, file ), sep = "\t" )


# tested genes of the differential expression analysis (columns: gene, symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, significance)
deg = as.data.frame( fread( paste0( degPath, "tables/DifferentialExpression/HOG_DEG_WithGeneSymbols.txt" ) ) )
pcGenes <- readLines( paste0( dataPath, "ReferenceFiles/GENCODE_v47/GENCODE_v47_ProteinCodingGenes.txt" ) )

# genes encoded in the mitochondrial genome
mtProteinGenes = c( "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6" )
mtAllGenes = c( mtProteinGenes, "MT-RNR1", "MT-RNR2", "MT-TA", "MT-TC", "MT-TD", "MT-TE", "MT-TF", "MT-TG", "MT-TH", "MT-TI", "MT-TK", "MT-TL1", "MT-TL2",
                "MT-TM", "MT-TN", "MT-TP", "MT-TQ", "MT-TR", "MT-TS1", "MT-TS2", "MT-TT", "MT-TV", "MT-TW", "MT-TY" )

# MitoCarta writes the mitochondrial genes without hyphen (MTND1), in our table they are called MT-ND1
unifySymbols <- function( x )
{
    x = toupper( as.character( x ) )
    oldNames = c( "MTND1", "MTND2", "MTND3", "MTND4", "MTND4L", "MTND5", "MTND6", "MTCO1", "MTCO2", "MTCO3", "MTCYB", "MTATP6", "MTATP8", "MTRNR1", "MTRNR2" )
    newNames = c( "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", "MT-ATP8", "MT-RNR1", "MT-RNR2" )
    cand = match( x, oldNames )
    x[ !is.na( cand ) ] = newNames[ cand[ !is.na( cand ) ] ]
    x
}

# MitoCarta 3.0 pathways (Rath et al., Nucleic Acids Res 2020; download see Data/ReferenceFiles/GeneSets/DOWNLOAD.txt)
# gmx format: one pathway per column, first row is the name, second row a description, below the genes
gmx = read.delim( paste0( dataPath, "ReferenceFiles/GeneSets/Human.MitoPathways3.0.gmx" ), header = FALSE, quote = "", comment.char = "", fill = TRUE, stringsAsFactors = FALSE )
mitoPathways = lapply( 1:ncol( gmx ), function( i ) {
  genes <- trimws( gmx[ -c( 1, 2 ), i ] )
  unique( genes[ !is.na( genes ) & genes != "" ] ) } )
names( mitoPathways ) <- trimws( as.character( gmx[ 1, ] ) )

targetPathway = "OXPHOS_subunits"


# running enrichment score of one gene set along the ranked gene list
computeRunningScore <- function( stats, pathwayGenes )
{
    genes = names( stats )
    hits = genes %in% pathwayGenes

    weightHit = abs( stats )
    weightHit[ !hits ] = 0
    running = cumsum( weightHit / sum( weightHit ) - ifelse( hits, 0, 1 / sum( !hits ) ) )

    # position of the maximal deviation from zero
    maxPos = which.max( running )
    minPos = which.min( running )
    peak = if( abs( running[ maxPos ] ) >= abs( running[ minPos ] ) ) maxPos else minPos

    # leading edge: genes of the set before the peak (or behind it for a negative score)
    position = 1:length( stats )
    leadingEdge = hits & if( running[ peak ] >= 0 ) position <= peak else position >= peak

    data.frame( rank = position, gene = genes, stat = as.numeric( stats ), hit = hits, running_es = running, leading_edge = leadingEdge )
}


# GSEA of all MitoCarta pathways, excludeMtGenes: remove the genes of the mitochondrial genome from ranking and pathways
runMitoGSEA <- function( excludeMtGenes, version )
{
    # ranking of protein-coding genes by the Wald statistic
    rankTab = deg[ deg$gene %in% pcGenes & !is.na( deg$padj ) & !is.na( deg$pvalue ) & !is.na( deg$stat ) & !is.na( deg$symbol ) & deg$symbol != "" & deg$stat != 0, ]
    rankTab$symbol = unifySymbols( rankTab$symbol )

    # ribosomal protein genes are not used
    isRibo = grepl( "^RP[LS]", rankTab$symbol )
    isMt = rankTab$symbol %in% mtAllGenes
    n_ribo = length( unique( rankTab$symbol[ isRibo ] ) )
    n_mt = if( excludeMtGenes ) length( unique( rankTab$symbol[ isMt ] ) ) else 0
    rankTab = rankTab[ !isRibo & !( excludeMtGenes & isMt ), ]

    # if a symbol occurs twice the entry with the larger statistic is taken
    rankTab = rankTab[ order( -abs( rankTab$stat ), rankTab$pvalue ), ]
    rankTab = rankTab[ !duplicated( rankTab$symbol ), ]
    stats = rankTab$stat
    names( stats ) = rankTab$symbol
    stats = sort( stats, decreasing = TRUE )

    # same gene names and filter for the pathways
    pathways = lapply( mitoPathways, function( genes ) {
      genes = unique( unifySymbols( genes ) )
      genes = genes[ !grepl( "^RP[LS]", genes ) ]
      if( excludeMtGenes ) genes = setdiff( genes, mtAllGenes )
      genes } )

    set.seed( 42 )
    fgseaRes = as.data.frame( fgseaMultilevel( pathways = pathways, stats = stats, minSize = 5, maxSize = 500, eps = 0 ) )
    fgseaRes$leadingEdge = sapply( fgseaRes$leadingEdge, paste, collapse = ";" )
    fgseaRes = fgseaRes[ order( fgseaRes$padj, -abs( fgseaRes$NES ) ), ]
    fgseaRes$version = version
    resRow = fgseaRes[ fgseaRes$pathway == targetPathway, ]

    # enrichment curve of the OXPHOS subunits
    targetGenes = pathways[[ targetPathway ]]
    overlap = intersect( targetGenes, names( stats ) )
    curve = computeRunningScore( stats, overlap )

    overview = data.frame( version = version,
                           pathway_reference_gene_count = length( targetGenes ),
                           ranked_overlap_gene_count = length( overlap ),
                           ranked_gene_count = length( stats ),
                           mt_pc_in_pathway_count = sum( mtProteinGenes %in% targetGenes ),
                           mt_pc_in_rank_count = sum( mtProteinGenes %in% names( stats ) ),
                           mt_pc_in_overlap_count = sum( mtProteinGenes %in% overlap ),
                           chr_m_symbols_removed_from_rank = n_mt,
                           rpl_rps_symbols_removed_from_rank = n_ribo,
                           ES = resRow$ES, NES = resRow$NES, pval = resRow$pval, padj = resRow$padj,
                           tested_mitopathway_count = nrow( fgseaRes ),
                           leading_edge_gene_count = sum( curve$leading_edge ) )

    # which genes of the original OXPHOS subunit list were used
    allTargetGenes = unique( unifySymbols( mitoPathways[[ targetPathway ]] ) )
    membership = data.frame( gene = allTargetGenes, version = version,
                             is_mt_pc = allTargetGenes %in% mtProteinGenes,
                             in_version_pathway = allTargetGenes %in% targetGenes,
                             in_rank = allTargetGenes %in% names( stats ),
                             in_ranked_pathway_overlap = allTargetGenes %in% overlap,
                             in_leading_edge = allTargetGenes %in% curve$gene[ curve$leading_edge ] )

    # statistics above 8 are cut in the plot
    rankPlot = data.frame( rank = 1:length( stats ), gene = names( stats ), stat = as.numeric( stats ), stat_plot = pmax( pmin( as.numeric( stats ), 8 ), -8 ) )

    list( overview = overview, membership = membership, fgseaRes = fgseaRes, curve = curve, rankPlot = rankPlot, resRow = resRow )
}


# Figure 2B: all protein-coding genes, Figure 2C: without the genes of the mitochondrial genome (do the nuclear encoded subunits alone show the same?)
resAll = runMitoGSEA( excludeMtGenes = FALSE, version = "mt-pc retained" )
resNuclear = runMitoGSEA( excludeMtGenes = TRUE, version = "mt-pc excluded" )

writeTable( rbind( resAll$overview, resNuclear$overview ), "GSEA_OXPHOS_Subunits_Overview.txt" )
writeTable( rbind( resAll$membership, resNuclear$membership ), "GSEA_OXPHOS_Subunits_GeneMembership.txt" )
writeTable( rbind( resAll$fgseaRes, resNuclear$fgseaRes ), "GSEA_MitoCarta_AllPathways.txt" )

# ranking and leading edge genes of the analysis with nuclear genes only
writeTable( resNuclear$rankPlot, "GSEA_NuclearGenes_WaldRanking.txt.gz" )
writeTable( resNuclear$curve[ resNuclear$curve$leading_edge, c( "gene", "rank", "stat" ) ], "GSEA_NuclearGenes_OXPHOS_LeadingEdge.txt" )

# data of the plots in 04_GSEA_MitoCarta_Plots.R
saveRDS( list( all = resAll[ c( "curve", "rankPlot", "resRow" ) ], nuclear = resNuclear[ c( "curve", "rankPlot", "resRow" ) ] ), file = paste0( objectPath, "GSEA_PlotData.rds" ) )
