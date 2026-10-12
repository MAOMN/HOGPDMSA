# change to your project directory
home_directory <- "path/to/project_directory/"

# results of the propeller test (made by 05_CellStateProportions_Analysis.R of Figure 3)
propPath = paste0( home_directory, "Figure_3/objects/" )

# output folder
outputPath = paste0( home_directory, "SupplementaryFigures/SupplementaryFigure_1A/" )
codePath = paste0( outputPath, "code/" )
figurePath = paste0( outputPath, "figures/" )

source( paste0( codePath, "saveFigure.R" ) )


# sizes of the figure, all in points (the files have the final size of the figure)
textSize = 10          # axis titles, panel titles
annotationSize = 8     # adjusted p-values above the lines
tickTextSize = 9       # numbers at the axes
letterSize = 13        # panel letters, bold
lineWidth = 0.5        # axes, ticks, lines above the groups
tickLength = 2
plotWidth = 91.25      # width and height of the plot area of one panel
plotHeight = 91
margins = c( 20, 39, 17, 6 )   # bottom, left, top, right of each panel
figureWidth = 4 * ( margins[ 2 ] + plotWidth + margins[ 4 ] )
cellHeight = plotHeight + margins[ 1 ] + margins[ 3 ]       # one panel with its margins
outerMargin = 4        # around the whole figure
# distances from the axis: axis title, numbers
titleDistance = 28
labelDistance = 4
# free space above the highest point for the two lines with the adjusted p-values (fraction of the range of the proportions)
headroom = 0.5
# position of the two lines above the highest point (fraction of the range of the proportions)
lineHeights = c( 0.12, 0.30 )
# distance between a line and its p-value, in points
valueDistance = 1.5
# point diameter of the figures of the propeller paper (13 x 6 inch, cex = 2), scaled to this figure
pointDiameter = 0.75 * 2 * 0.66 * 12 * 545 / ( 13 * 72 )   # fixed, the points do not change with the size of the figure

# labels of the y-axis: no trailing zeros after the decimal point, zero is written as 0
yLabels <- function( x )
{
    res = format( x, trim = TRUE, scientific = FALSE )
    res = sub( "[.]0+$", "", res )
    res = sub( "([.][0-9]*[1-9])0+$", "\\1", res )
    res[ abs( x ) < .Machine$double.eps^0.5 ] = "0"
    res
}

states = c( "OPC", paste0( "OL", 1:6 ) )
groups = c( "HC", "MSA", "PD" )
# same colors as in the UMAP of Figure 3A
diagnosis_color = c( HC = "#D0CDD6", MSA = "#E06666", PD = "#6FA8DC" )

# results of the two propeller tests (each disease against HC, adjusted for age)
# list with: props (output of getTransformedProps), donor_info, result (output of propeller.ttest)
prop = list( MSA = readRDS( paste0( propPath, "Propeller_MSA_vs_HC.rds" ) ), PD = readRDS( paste0( propPath, "Propeller_PD_vs_HC.rds" ) ) )

# proportions of every donor: the HC donors are the same in both tests, they are taken once
donorTab = c()
for( disease in c( "MSA", "PD" ) )
{
    info = prop[[ disease ]]$donor_info
    for( state in states )
    {
        y = as.numeric( prop[[ disease ]]$props$Proportions[ state, as.character( info$sample_id ) ] )
        donorTab = rbind( donorTab, data.frame( state = state, donor = as.character( info$sample_id ), group = as.character( info$Diagnosis ), proportion = y ) )
    }
}
donorTab = donorTab[ !duplicated( donorTab[ , c( "state", "donor" ) ] ), ]
donorTab$group = factor( donorTab$group, levels = groups )

# adjusted p-value of the test, written with three decimals and without trailing zeros
pLabel <- function( p ) as.character( round( p, 3 ) )

#
# One panel per cell state with the proportions of all donors: HC, MSA and PD next to each other
# (stripchart and title as in the result figures of the propeller paper, Phipson et al. 2022, RealDataAnalysis.Rmd)
# the adjusted p-value of each test is written above a line over the two compared groups
#
drawPanels <- function( letter )
{
    for( state in states )
    {
        d = donorTab[ donorTab$state == state, ]
        y = d$proportion
        r = diff( range( y ) )
        stripchart( proportion ~ group, data = d, vertical = TRUE, pch = 16, method = "jitter",
                    col = diagnosis_color[ groups ], cex = pointDiameter / ( 0.75 * textSize ), ylab = "Proportions",
                    ylim = c( min( y ), max( y ) + headroom * r ), cex.lab = 1, axes = FALSE, frame.plot = FALSE )
        # axes at the left and at the bottom only, no frame at the top and at the right
        box( which = "plot", bty = "l", lwd = lineWidth * 96 / 72 )
        axis( 1, at = 1:3, labels = groups, cex.axis = tickTextSize / textSize, lwd = 0, lwd.ticks = lineWidth * 96 / 72 )
        # numbers without trailing zeros after the decimal point (0.1 instead of 0.10), zero is written as 0
        ticks = axTicks( 2 )
        axis( 2, at = ticks, labels = yLabels( ticks ), cex.axis = tickTextSize / textSize, lwd = 0, lwd.ticks = lineWidth * 96 / 72, las = 1 )
        title( state, cex.main = 1, font.main = 1, adj = 0, line = 4 / linePt )

        # line over HC and MSA (lower) and over HC and PD (upper), the adjusted p-value of the test above it
        unitPerPoint = ( 1 + headroom ) * r / plotHeight
        lines = list( list( to = 2, disease = "MSA", height = lineHeights[ 1 ] ), list( to = 3, disease = "PD", height = lineHeights[ 2 ] ) )
        for( l in lines )
        {
            yLine = max( y ) + l$height * r
            segments( 1, yLine, l$to, yLine, lwd = lineWidth * 96 / 72, lend = 1, xpd = NA )
            text( ( 1 + l$to ) / 2, yLine + valueDistance * unitPerPoint, pLabel( prop[[ l$disease ]]$result[ state, "FDR" ] ),
                  adj = c( 0.5, 0 ), cex = annotationSize / textSize, xpd = NA )
        }

        # panel letter at the upper left corner of the figure (not for the single panel)
        if( state == states[ 1 ] && !is.null( letter ) )
            text( grconvertX( 0, "nic", "user" ), grconvertY( 1, "nfc", "user" ), letter, adj = c( 0, 1 ), font = 2, cex = letterSize / textSize, xpd = NA )
    }
}

setSizes <- function()
{
    par( family = "Arial", ps = textSize, cex = 1, omi = rep( outerMargin, 4 ) / 72, mai = margins / 72, col.axis = "black", col.lab = "black" )
    # height of one margin line in points, to give the distances in points
    linePt <<- 72 * par( "mai" )[ 1 ] / par( "mar" )[ 1 ]
    par( mgp = c( titleDistance, labelDistance, 0 ) / linePt, tcl = -tickLength / linePt )
}

# Supplementary Figure 1A: seven cell states in two rows
drawFigure <- function( letter )
{
    layout( matrix( c( 1:4, 5:7, 0 ), ncol = 4, byrow = TRUE ) )
    setSizes()
    drawPanels( letter )
}

# figure with the panel letter and the panel without it
saveFigure( function() drawFigure( "A" ), paste0( figurePath, "SupplementaryFigure_1A" ), width = ( figureWidth + 2 * outerMargin ) / 72, height = ( 2 * cellHeight + 2 * outerMargin ) / 72, dpi = 600 )
saveFigure( function() drawFigure( NULL ), paste0( figurePath, "Panels/SupplementaryFigure_1A" ), width = ( figureWidth + 2 * outerMargin ) / 72, height = ( 2 * cellHeight + 2 * outerMargin ) / 72, dpi = 600 )
