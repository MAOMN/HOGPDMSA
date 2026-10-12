library(tradeSeq)
library(slingshot)
library(SingleCellExperiment)
library(Matrix)
library(edgeR)
library(BiocParallel)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# trajectories and counts of the genes used for the models (Figure 4)
trajPath = paste0( home_directory, "Figure_4/" )

# output folder
outputPath = paste0( home_directory, "Figure_5/" )
tablePath = paste0( outputPath, "tables/TwoLineages/" )
objectPath = paste0( outputPath, "objects/" )


# comparison of each disease with HC along the trajectories
# - separate models for the MSA + HC cells and for the PD + HC cells
# - lineage 1 (to OL5) and lineage 3 (to OL3), the OL6 cluster is not part of this comparison

# slingshot result (made by 01_Trajectories_Slingshot_Analysis.R of Figure 4)
# list with: slingshot_sce, lineages, pseudotime, cell_weights, cell_metadata
trajectory <- readRDS( paste0( trajPath, "objects/OPC_OL_Slingshot.rds" ) )
cells <- trajectory$cell_metadata

# counts of the genes that are modelled, genes x cells (made by 07_fitGAM_Analysis.R of Figure 4)
counts <- readRDS( paste0( trajPath, "objects/fitGAM_Counts.rds" ) )

lineagesUsed = c( 1, 3 )
cellWeights <- trajectory$cell_weights[ , lineagesUsed ]
pseudotime <- slingPseudotime( trajectory$slingshot_sce, na = FALSE )[ , lineagesUsed ]
colnames( cellWeights ) <- colnames( pseudotime ) <- c( "L1", "L3" )

# cells without OL6 that have a weight on at least one of the two lineages
keep <- cells$final_cluster_label != "OL6" & rowSums( cellWeights ) > 0
cells <- cells[ keep, ]
counts <- counts[ , keep ]
cellWeights <- cellWeights[ keep, ]
pseudotime <- pseudotime[ keep, ]

# two genes have no counts at all in the HC and MSA cells, they are removed (also for PD, so that both comparisons use the same genes)
noCounts <- rowSums( counts[ , as.character( cells$diagnosis ) %in% c( "HC", "MSA" ) ] ) == 0
genes = rownames( counts )[ !noCounts ]

modelDir = c( MSA = "fitGAM_MSA_vs_HC_TwoLineages", PD = "fitGAM_PD_vs_HC_TwoLineages" )
seeds = c( MSA = 202608311, PD = 202608312 )
nKnots = 6

rangeTab = c()
for( disease in c( "MSA", "PD" ) )
{
    selected = as.character( cells$diagnosis ) %in% c( "HC", disease )
    countsPair = counts[ genes, selected ]
    pt = pseudotime[ selected, ]
    wt = cellWeights[ selected, ]
    cellsPair = cells[ selected, ]
    conditions = factor( as.character( cellsPair$diagnosis ), levels = c( "HC", disease ) )

    # offset of each cell: library size with TMM normalization factor
    libSize = colSums( countsPair )
    normFactors = calcNormFactors( countsPair, lib.size = libSize, method = "TMM" )
    offset = log( libSize * normFactors )
    names( offset ) = colnames( countsPair )

    # the models are fitted in batches of 100 genes
    batchTab = data.frame( gene = genes, gene_rank = 1:length( genes ), detected_cells = rowSums( countsPair > 0 ), total_counts = rowSums( countsPair ),
                           batch_id = ceiling( 1:length( genes ) / 100 ) )
    fwrite( batchTab, file = paste0( tablePath, "fitGAM_GeneBatches_", disease, "_vs_HC.txt.gz" ), sep = "\t" )

    saveRDS( list( counts = countsPair, pseudotime = pt, cell_weights = wt, cell_metadata = cellsPair, conditions = conditions, offset = offset ),
             file = paste0( objectPath, "fitGAM_Inputs_", disease, "_vs_HC.rds" ) )

    # fitGAM needs several days for all genes, the models are saved in batches of 100 genes
    modelPath = paste0( objectPath, modelDir[ disease ], "/" )
    dir.create( modelPath, showWarnings = FALSE )
    for( b in unique( batchTab$batch_id ) )
    {
        g <- batchTab$gene[ batchTab$batch_id == b ]
        set.seed( seeds[ disease ] )
        fit <- fitGAM( counts = countsPair[ g, ], pseudotime = pt, cellWeights = wt, conditions = conditions, offset = offset, nknots = nKnots,
                       parallel = TRUE, BPPARAM = MulticoreParam( 10, RNGseed = seeds[ disease ] ), control = mgcv::gam.control( maxit = 100 ), sce = TRUE )
        saveRDS( fit, file = paste0( modelPath, "fit_sce_batch_", sprintf( "%04d", b ), ".rds" ) )
    }

    #
    # Part of the lineages that is compared between disease and HC
    #
    # the comparison is restricted to the pseudotime range of the oligodendrocytes (OPCs are at the start of both lineages)
    # the range is given by the knots of the model: first knot inside the range of the OL cells until the last knot inside
    # the knots are the same in all batches
    knots = as.numeric( metadata( readRDS( paste0( modelPath, "fit_sce_batch_0001.rds" ) ) )$tradeSeq$knots )
    isOL = cellsPair$final_cluster_label != "OPC"
    for( i in 1:2 )
    {
        olCells = isOL & is.finite( pt[ , i ] ) & wt[ , i ] > 0
        ptMin = min( pt[ olCells, i ] )
        ptMax = max( pt[ olCells, i ] )
        knotStart = which( knots >= ptMin )[ 1 ]
        knotEnd = tail( which( knots <= ptMax ), 1 )
        rangeTab = rbind( rangeTab, data.frame( disease = disease, model_lineage_index = i, original_lineage_index = lineagesUsed[ i ], lineage_id = c( "L1", "L3" )[ i ], endpoint = c( "OL5", "OL3" )[ i ],
                                                positive_weight_OL_cells = sum( olCells ), ol_pseudotime_min = ptMin, ol_pseudotime_max = ptMax,
                                                test_knot_start = knotStart, test_knot_end = knotEnd, test_time_start = knots[ knotStart ], test_time_end = knots[ knotEnd ],
                                                actual_knots = paste( knots, collapse = ";" ) ) )
    }
}

fwrite( rangeTab, file = paste0( tablePath, "Disease_Dynamics_ComparedPseudotimeRange.txt" ), sep = "\t" )
saveRDS( rangeTab, file = paste0( objectPath, "Disease_Dynamics_ComparedPseudotimeRange.rds" ) )
