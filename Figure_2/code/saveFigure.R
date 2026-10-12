library(ragg)
library(svglite)
library(xml2)

# every figure is saved three times: pdf, png and svg
# the svg files are the ones used to put the panels together, the text in them stays editable
# font is always Arial

# without this R opens Rplots.pdf as soon as a text size is measured and no device is open
options( device = function( ... ) ragg::agg_capture( width = 1500, height = 1500, res = 300 ) )
grDevices::pdfFonts( Arial = grDevices::pdfFonts( "ArialMT" )[[ 1 ]] )

# pdf: quartz on macOS (the letters are placed at the exact widths of the font, every text stays one piece in Illustrator),
# cairo_pdf on other systems
openPdf <- function( file, width, height )
{
    if( capabilities( "aqua" ) ) quartz( file = file, type = "pdf", width = width, height = height, family = "Arial", bg = "white" )
    else cairo_pdf( file, width = width, height = height, family = "Arial", bg = "white" )
}

#
# Save figure as pdf, png and svg
#
# p - ggplot object or a function that draws the figure
# fileName - output file without ending
# pngFromPdf - make the png out of the pdf (for figures with very many points)
saveFigure <- function( p, fileName, width, height, dpi = 300, pngFromPdf = FALSE )
{
    drawFigure <- function()
    {
        #positions of repelled labels are random, keep them the same in all three files
        set.seed( 42 )
        if( is.function( p ) ) p() else print( p )
    }

    openPdf( paste0( fileName, ".pdf" ), width, height )
    drawFigure()
    dev.off()

    if( pngFromPdf )
    {
        system2( "pdftocairo", c( "-png", "-singlefile", "-r", dpi, shQuote( paste0( fileName, ".pdf" ) ), shQuote( fileName ) ) )
    }
    else
    {
        agg_png( paste0( fileName, ".png" ), width = width, height = height, units = "in", res = dpi, background = "white" )
        drawFigure()
        dev.off()
    }

    svgFile = paste0( fileName, ".svg" )
    svglite( svgFile, width = width, height = height, bg = "white", standalone = TRUE, always_valid = TRUE, fix_text_size = FALSE,
             system_fonts = list( sans = "Arial", serif = "Arial", mono = "Arial" ) )
    drawFigure()
    dev.off()

    #clipping paths and masks make the svg hard to edit afterwards, remove them
    svg = read_xml( svgFile )
    xml_remove( xml_find_all( svg, "//*[local-name()='clipPath' or local-name()='mask']" ) )
    for( a in c( "clip-path", "mask" ) )
    {
        nodes = xml_find_all( svg, paste0( "//*[@", a, "]" ) )
        xml_attr( nodes, a ) <- NULL
    }
    #svglite sets the minus sign and some spaces of expressions in the font Symbol, all text of the figures is Arial
    styleNodes = xml_find_all( svg, "//*[contains(@style, 'Symbol')]" )
    if( length( styleNodes ) > 0 ) xml_attr( styleNodes, "style" ) <- gsub( "font-family: *\"Symbol\"", "font-family: \"Arial\"", xml_attr( styleNodes, "style" ) )
    write_xml( svg, svgFile )
}
