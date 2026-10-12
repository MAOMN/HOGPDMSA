library(tximport)
library(DESeq2)
library(org.Hs.eg.db)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_1/" )
tablePath = paste0( outputPath, "tables/" )
objectPath = paste0( outputPath, "objects/" )

writeTable <- function( x, file, quote = FALSE ) write.table( x, file = paste0( tablePath, file ), col.names = TRUE, row.names = FALSE, quote = quote, sep = "\t" )


# overview of the HOG samples (columns: sample, condition, replicate)
sample_info = read.delim( paste0( dataPath, "HOG_RNAseq/HOG_SampleOverview.txt" ) )
sample_info$condition <- factor( sample_info$condition, levels = c( "Control", "asyn_PFF" ) )
rownames( sample_info ) <- sample_info$sample

# Salmon output (one folder per sample) and the assignment of transcripts to genes of GENCODE v47 (gene ids without version number)
quantFiles = setNames( paste0( dataPath, "HOG_RNAseq/HOG_SalmonQuant/", sample_info$sample, "/quant.sf" ), sample_info$sample )
tx2gene = read.delim( paste0( dataPath, "ReferenceFiles/GENCODE_v47/GENCODE_v47_tx2gene.txt" ) )
tx2gene$GENEID <- sub( "[.][0-9]+$", "", tx2gene$GENEID )

# transcript counts summarized to genes, genes with almost no reads are removed
txi <- tximport( quantFiles, type = "salmon", tx2gene = tx2gene )
dds <- DESeqDataSetFromTximport( txi, colData = sample_info, design = ~ condition )
dds <- dds[ rowSums( counts( dds ) ) >= 10, ]

# differential expression alpha-syn PFF vs control, rlog values for PCA and sample distances
dds <- DESeq( dds, fitType = "local" )
res <- results( dds, contrast = c( "condition", "asyn_PFF", "Control" ), alpha = 0.05 )
rld <- rlog( dds, blind = FALSE )

saveRDS( dds, file = paste0( objectPath, "HOG_DESeq2_dds.rds" ) )
saveRDS( res, file = paste0( objectPath, "HOG_DESeq2_results.rds" ) )
saveRDS( rld, file = paste0( objectPath, "HOG_rlog.rds" ) )

resTable = data.frame( gene = rownames( res ), as.data.frame( res ), row.names = NULL )
writeTable( resTable, "DifferentialExpression/HOG_DESeq2_AllGenes.txt" )

# normalization factors (tximport gives one factor per gene and sample)
nf <- normalizationFactors( dds )
writeTable( data.frame( gene = rownames( nf ), nf ), "DifferentialExpression/HOG_DESeq2_NormalizationFactors.txt" )
writeTable( data.frame( sample = colnames( dds ), minimum = apply( nf, 2, min ), median = apply( nf, 2, median ), maximum = apply( nf, 2, max ),
                        normalization = "gene-specific tximport offset" ), "DifferentialExpression/HOG_DESeq2_NormalizationFactors_Summary.txt" )

#
### Quality control: PCA on the 500 genes with highest variance and overview of the analysis
#

rlogValues <- assay( rld )
topVar <- order( rowVars( rlogValues ), decreasing = TRUE )[ 1:500 ]
pca <- prcomp( t( rlogValues[ topVar, ] ) )
varExplained <- pca$sdev^2 / sum( pca$sdev^2 )
df_pca = data.frame( PC1 = pca$x[ , 1 ], PC2 = pca$x[ , 2 ], sample = colnames( rld ), cond = sample_info$condition,
                     PC1_var_pct = 100 * varExplained[ 1 ], PC2_var_pct = 100 * varExplained[ 2 ] )
writeTable( df_pca, "DifferentialExpression/HOG_PCA_Coordinates.txt" )

libSize <- colSums( counts( dds ) )
sigRes <- res[ !is.na( res$padj ) & res$padj < 0.05, ]
writeTable( data.frame( dataset = "HOG",
                        n_samples = ncol( dds ),
                        n_groups = length( unique( sample_info$condition ) ),
                        design = "~condition",
                        contrast = "condition; asyn_PFF; Control",
                        n_genes_raw = nrow( dds ),
                        n_genes_tested = sum( !is.na( res$padj ) ),
                        n_genes_NA = sum( is.na( res$padj ) ),
                        n_DEGs_padj05 = nrow( sigRes ),
                        n_Up = sum( sigRes$log2FoldChange > 0 ),
                        n_Down = sum( sigRes$log2FoldChange < 0 ),
                        median_libsize = median( libSize ),
                        min_libsize = min( libSize ),
                        max_libsize = max( libSize ),
                        pc1_var_pct = round( 100 * varExplained[ 1 ], 2 ),
                        pc2_var_pct = round( 100 * varExplained[ 2 ], 2 ),
                        rlog_transform = "rlog (N<30)" ), "DifferentialExpression/HOG_QC_Summary.txt" )

#
### Table of the tested genes with gene symbols
#

resTable$symbol <- mapIds( org.Hs.eg.db, keys = resTable$gene, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first" )[ resTable$gene ]

# genes without p-value were not tested (no reads in one group or outlier)
resTable <- resTable[ !is.na( resTable$pvalue ), ]
resTable$significance <- factor( ifelse( is.na( resTable$padj ) | resTable$padj >= 0.05, "NS", ifelse( resTable$log2FoldChange > 0, "Up", "Down" ) ),
                                 levels = c( "Up", "Down", "NS" ) )
resTable <- resTable[ , c( "gene", "symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significance" ) ]
writeTable( resTable[ order( resTable$padj, na.last = TRUE ), ], "DifferentialExpression/HOG_DEG_WithGeneSymbols.txt" )

#
### Data of the volcano plot (Figure 1C): protein-coding genes only, genes without adjusted p-value are shown at the bottom
#

pcGenes <- readLines( paste0( dataPath, "ReferenceFiles/GENCODE_v47/GENCODE_v47_ProteinCodingGenes.txt" ) )
volcanoData <- resTable[ resTable$gene %in% pcGenes, ]
volcanoData$logFDR <- -log10( pmax( volcanoData$padj, 1e-300 ) )
volcanoData$logFDR[ is.na( volcanoData$logFDR ) ] <- 0

# genes with a label: the 5 most significant ones and the three genes named in the Results text
# labels of genes with a very large fold change are left out
sigGenes <- volcanoData[ volcanoData$significance != "NS", ]
sigGenes <- sigGenes[ order( sigGenes$padj ), ]
lfcLimit <- max( 3, min( 8, ceiling( quantile( abs( sigGenes$log2FoldChange ), 0.95, names = FALSE ) + 0.3 ) ) )
topGenes <- sigGenes$symbol[ abs( sigGenes$log2FoldChange ) <= lfcLimit ][ 1:5 ]
labelGenes <- volcanoData[ match( unique( c( topGenes, "CXCL8", "IL1A", "IL1B" ) ), volcanoData$symbol ), ]
stopifnot( !anyNA( labelGenes$gene ), all( labelGenes$significance != "NS" ) )
writeTable( labelGenes[ , 1:9 ], "DifferentialExpression/HOG_VolcanoPlot_LabelledGenes.txt" )

# data of the plots in 02_DifferentialExpression_Plots.R
saveRDS( list( pca = df_pca, volcano = volcanoData, labels = labelGenes ), file = paste0( objectPath, "HOG_PlotData.rds" ) )
