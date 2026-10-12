# change to your project directory
home_directory <- "path/to/project_directory/"

# data folder
dataPath = paste0( home_directory, "Data/" )

# output folder
outputPath = paste0( home_directory, "Figure_4/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )

source( paste0( codePath, "saveFigure.R" ) )


# the plot function of tradeSeq (tradeSeq:::plot_evalutateK_results), only the numbers at the y axes are written without a zero at the end (0.9990 is written 0.999)
axisNumbers <- function( v ) format( v, scientific = FALSE, trim = TRUE, drop0trailing = TRUE )
# all ticks, the numbers at every "every"-th tick (the default function leaves out numbers that would overlap)
axisY <- function( every = 1 )
{
    ticks = graphics::axTicks( 2 )
    graphics::axis( 2, at = ticks, labels = FALSE )
    shown = ticks[ seq( 1, length( ticks ), by = every ) ]
    graphics::axis( 2, at = shown, labels = axisNumbers( shown ), tick = FALSE )
}
plotEvaluateK <- function( aicMat, k, aicDiff = 2 )
{
    init_shape <- graphics::par()$mfrow
    graphics::par( mfrow = c( 1, 4 ) )
    devs <- matrix( NA, nrow = nrow( aicMat ), ncol = length( k ) )
    for( ii in seq_len( nrow( aicMat ) ) ) devs[ ii, ] <- aicMat[ ii, ] - mean( aicMat[ ii, ] )
    graphics::boxplot( devs, ylab = "Deviation from genewise average AIC", xlab = "Number of knots", xaxt = "n", yaxt = "n" )
    graphics::axis( 1, at = seq_len( length( k ) ), labels = k )
    axisY()
    plot( x = k, y = colMeans( aicMat, na.rm = TRUE ), type = "b", ylab = "Average AIC", xlab = "Number of knots", yaxt = "n" )
    axisY()
    plot( x = k, y = colMeans( aicMat / aicMat[ , 1 ], na.rm = TRUE ), type = "b", ylab = "Relative AIC", xlab = "Number of knots", yaxt = "n" )
    axisY( every = 2 )
    aicRange <- apply( apply( aicMat, 1, range ), 2, diff )
    varID <- which( aicRange > aicDiff )
    if( length( varID ) > 0 )
    {
        aicMatSub <- aicMat[ varID, ]
        tab <- table( k[ apply( aicMatSub, 1, which.min ) ] )
        graphics::barplot( tab, xlab = "Number of knots", ylab = "# Genes with optimal k" )
    }
    graphics::par( mfrow = init_shape )
}


#
### Diagnostic plots of evaluateK (the plot function of tradeSeq on the saved AIC values, made by 05_evaluateK_Analysis.R):
### deviation from the mean AIC of the gene, mean AIC, AIC relative to 3 knots, number of genes with the lowest AIC at each number of knots (genes with an AIC range above 2)
#

for( i in 1:3 )
{
    aic = readRDS( paste0( dataPath, "SavedModels/evaluateK/L", i, "_k3_10.rds" ) )$aic
    k = as.integer( sub( "k: ", "", colnames( aic ) ) )
    saveFigure( function() plotEvaluateK( aicMat = aic, k = k, aicDiff = 2 ), paste0( figurePath, "fitGAM_NumberOfKnots_evaluateK_L", i ), width = 11, height = 3.2 )
}
