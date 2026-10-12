library(DESeq2)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
degPath = paste0( home_directory, "Figure_1/" )
outputPath = paste0( home_directory, "Figure_2/" )
tablePath = paste0( outputPath, "tables/MitochondrialGenes/" )
objectPath = paste0( outputPath, "objects/" )

writeTable <- function( x, file ) fwrite( x, file = paste0( tablePath, file ), sep = "\t" )


# the 13 protein-coding genes of the mitochondrial genome (columns: symbol, ensg)
mtGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GENCODE_v47/GENCODE_v47_MitochondrialProteinCodingGenes.txt" ) ) )

# rlog expression values and test results of the HOG cells (01_DifferentialExpression_Analysis.R of Figure 1)
rld <- readRDS( paste0( degPath, "objects/HOG_rlog.rds" ) )
res <- as.data.frame( readRDS( paste0( degPath, "objects/HOG_DESeq2_results.rds" ) ) )

# controls first, then alpha-syn treated samples
sample_info <- as.data.frame( colData( rld ) )
samples <- c( rownames( sample_info )[ sample_info$condition == "Control" ], rownames( sample_info )[ sample_info$condition == "asyn_PFF" ] )
group <- sapply( samples, function( s ) if ( sample_info[ s, "condition" ] == "Control" ) "Control" else "alpha-syn" )
group <- factor( group, levels = c( "Control", "alpha-syn" ) )

expr <- assay( rld )[ mtGenes$ensg, samples ]
rownames( expr ) <- mtGenes$symbol

# z-score of each gene over the six samples
z_expr <- ( expr - rowMeans( expr ) ) / apply( expr, 1, sd )

# test results of the 13 genes
mtRes <- res[ mtGenes$ensg, ]
mtRes$star <- sapply( mtRes$padj, function( x ) if ( is.na( x ) ) "" else if ( x < 0.001 ) "***" else if ( x < 0.01 ) "**" else if ( x < 0.05 ) "*" else "" )
mtRes <- cbind( mtGenes[ , c( "ensg", "symbol" ) ], mtRes )
writeTable( mtRes, "MitochondrialGenes_DESeq2_Statistics.txt" )

# z-scores as long table (columns: symbol, sample, condition, z)
writeTable( data.frame( symbol = rep( rownames( z_expr ), ncol( z_expr ) ),
                        sample = rep( colnames( z_expr ), each = nrow( z_expr ) ),
                        condition = rep( as.character( group ), each = nrow( z_expr ) ),
                        z = as.vector( z_expr ) ), "MitochondrialGenes_ZScores.txt" )

# data of the heatmap in 02_MitochondrialGenes_Plots.R
saveRDS( list( z_expr = z_expr, group = group, mtRes = mtRes ), file = paste0( objectPath, "MitochondrialGenes_PlotData.rds" ) )
