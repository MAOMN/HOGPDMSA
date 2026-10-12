library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# results of the clustering (Figure 3)
clusterPath = paste0( home_directory, "Figure_3/" )

# output folder
outputPath = paste0( home_directory, "Figure_4/" )
tablePath = paste0( outputPath, "tables/Trajectories/" )
objectPath = paste0( outputPath, "objects/" )


# table of the clustered OPCs and oligodendrocytes (made by 03_OPC_OL_Clustering_Analysis.R of Figure 3), only cells of the kept clusters
cells = as.data.frame( fread( paste0( clusterPath, "tables/OPC_OL_Clustering/OPC_OL_CellTable.txt.gz" ) ) )
cells <- cells[ cells$final_retained, c( "cell_id", "sample_id", "diagnosis", "final_cluster_label" ) ]
cells$diagnosis <- factor( cells$diagnosis, levels = c( "HC", "MSA", "PD" ) )
rownames( cells ) <- cells$cell_id

# harmony corrected principal components (made by 03_OPC_OL_Clustering_Analysis.R of Figure 3)
harmony <- readRDS( paste0( clusterPath, "objects/OPC_OL_HarmonyEmbedding.rds" ) )
harmony <- harmony[ cells$cell_id, ]


#
### Trajectories with slingshot: in the space of the first three components, start in the OPC cluster
#

pc3 <- harmony[ , 1:3 ]
colnames( pc3 ) <- paste0( "HarmonyPC", 1:3 )

# slingshot only needs the coordinates and the clusters, the expression matrix is left empty
sce <- SingleCellExperiment( assays = list( dummy = matrix( 0, 1, nrow( cells ), dimnames = list( "dummy", cells$cell_id ) ) ), colData = DataFrame( cells ) )
reducedDim( sce, "HARMONY_PCA3" ) <- pc3

sce <- slingshot( sce, clusterLabels = "final_cluster_label", reducedDim = "HARMONY_PCA3", start.clus = "OPC", extend = "n", stretch = 0, approx_points = 150 )

lineages <- slingLineages( sce )

pseudotime <- slingPseudotime( sce )
cellWeights <- slingCurveWeights( sce )
colnames( pseudotime ) <- colnames( cellWeights ) <- paste0( "Lineage_", 1:length( lineages ) )

saveRDS( list( slingshot_sce = sce, lineages = lineages, pseudotime = pseudotime, cell_weights = cellWeights, cell_metadata = cells ),
         file = paste0( objectPath, "OPC_OL_Slingshot.rds" ) )

# clusters along each lineage
lineageTab = c()
for( i in 1:length( lineages ) )
{
    lineageTab = rbind( lineageTab, data.frame( lineage = paste( "Lineage", i ), path_order = 1:length( lineages[[ i ]] ), cluster = lineages[[ i ]],
                                                root_cluster = lineages[[ i ]][ 1 ], terminal_cluster = tail( lineages[[ i ]], 1 ) ) )
}
fwrite( lineageTab, file = paste0( tablePath, "Trajectories_LineagePaths.txt" ), sep = "\t" )


#
### Tables for the plots of the trajectories
#

# each cell is shown on the lineage where it has the largest weight
assigned = max.col( cellWeights, ties.method = "first" )
cellTab = cells
cellTab$assigned_lineage = paste( "Lineage", assigned )
cellTab$assigned_lineage_end_cluster = sapply( lineages, tail, 1 )[ assigned ]
cellTab$assigned_weight = cellWeights[ cbind( 1:nrow( cells ), assigned ) ]
cellTab$assigned_pseudotime = pseudotime[ cbind( 1:nrow( cells ), assigned ) ]
cellTab$PCA3_1 = pc3[ , 1 ]
cellTab$PCA3_2 = pc3[ , 2 ]
cellTab$PCA3_3 = pc3[ , 3 ]

# points of the fitted curves, ordered along the curve
curveTab = c()
for( i in 1:length( lineages ) )
{
    curve = slingCurves( sce )[[ i ]]
    xyz = as.data.frame( curve$s[ curve$ord, ] )
    colnames( xyz ) = c( "PCA3_1", "PCA3_2", "PCA3_3" )

    #remove points that are identical to the point before
    xyz$segment_length = c( 0, sqrt( diff( xyz$PCA3_1 )^2 + diff( xyz$PCA3_2 )^2 + diff( xyz$PCA3_3 )^2 ) )
    xyz = xyz[ c( TRUE, xyz$segment_length[ -1 ] > 1e-12 ), ]

    #position on the curve in units of pseudotime
    xyz$curve_arc = cumsum( c( 0, sqrt( diff( xyz$PCA3_1 )^2 + diff( xyz$PCA3_2 )^2 + diff( xyz$PCA3_3 )^2 ) ) )
    lambda = curve$lambda[ is.finite( curve$lambda ) ]
    xyz$curve_lambda = min( lambda ) + xyz$curve_arc / max( xyz$curve_arc ) * ( max( lambda ) - min( lambda ) )

    xyz$lineage = paste0( "Lineage ", i )
    xyz$terminal_cluster = tail( lineages[[ i ]], 1 )
    xyz$curve_order = 1:nrow( xyz )
    xyz$lambda_min = min( lambda )
    xyz$lambda_max = max( lambda )
    curveTab = rbind( curveTab, xyz[ , c( "lineage", "terminal_cluster", "curve_order", "curve_lambda", "PCA3_1", "PCA3_2", "PCA3_3", "curve_arc", "lambda_min", "lambda_max", "segment_length" ) ] )
}

# position of each cell on its curve (coordinates of the curve at the pseudotime of the cell)
attachedTab = cellTab[ , c( "cell_id", "sample_id", "diagnosis", "final_cluster_label", "assigned_lineage", "assigned_lineage_end_cluster", "assigned_weight", "assigned_pseudotime" ) ]
attachedTab$attached_PCA3_1 = NA
attachedTab$attached_PCA3_2 = NA
attachedTab$attached_PCA3_3 = NA
for( lin in unique( curveTab$lineage ) )
{
    cand = which( attachedTab$assigned_lineage == lin )
    curve = curveTab[ curveTab$lineage == lin, ]
    curve = curve[ order( curve$curve_lambda ), ]
    lambda = pmax( min( curve$curve_lambda ), pmin( max( curve$curve_lambda ), attachedTab$assigned_pseudotime[ cand ] ) )
    for( k in 1:3 )
    {
        attachedTab[ cand, paste0( "attached_PCA3_", k ) ] = approx( curve$curve_lambda, curve[ , paste0( "PCA3_", k ) ], xout = lambda, rule = 2, ties = "ordered" )$y
    }
}

fwrite( cellTab, file = paste0( tablePath, "Trajectories_Cells.txt.gz" ), sep = "\t" )
fwrite( curveTab, file = paste0( tablePath, "Trajectories_Curves.txt" ), sep = "\t" )
fwrite( attachedTab, file = paste0( tablePath, "Trajectories_CellsOnCurves.txt.gz" ), sep = "\t" )


#
### Three-dimensional UMAP to show the trajectories
#

# same input as for the two-dimensional UMAP (harmony components 1 to 30), cells sorted by their id
ids = sort( cells$cell_id )
umap3d <- Embeddings( RunUMAP( harmony[ ids, 1:30 ], n.components = 3 ) )
umap3d <- umap3d[ cells$cell_id, ]

umapTab = data.frame( cells, UMAP3D_1 = umap3d[ , 1 ], UMAP3D_2 = umap3d[ , 2 ], UMAP3D_3 = umap3d[ , 3 ] )
fwrite( umapTab, file = paste0( tablePath, "Trajectories_UMAP3D_Cells.txt.gz" ), sep = "\t" )

# read the saved coordinates again and map the slingshot curves into the UMAP
umapTab = as.data.frame( fread( paste0( tablePath, "Trajectories_UMAP3D_Cells.txt.gz" ) ) )
xyz <- as.matrix( umapTab[ , c( "UMAP3D_1", "UMAP3D_2", "UMAP3D_3" ) ] )
rownames( xyz ) <- umapTab$cell_id

embedded <- embedCurves( sce, xyz, approx_points = 250 )

umapCurveTab = c()
for( i in 1:length( lineages ) )
{
    curve = slingCurves( embedded )[[ i ]]
    pts = as.data.frame( curve$s[ curve$ord, ] )
    colnames( pts ) = c( "UMAP3D_1", "UMAP3D_2", "UMAP3D_3" )
    pts$segment_length = c( 0, sqrt( diff( pts$UMAP3D_1 )^2 + diff( pts$UMAP3D_2 )^2 + diff( pts$UMAP3D_3 )^2 ) )
    pts = pts[ c( TRUE, pts$segment_length[ -1 ] > 1e-12 ), ]
    pts$lineage = paste( "Lineage", i )
    pts$curve_order = 1:nrow( pts )
    umapCurveTab = rbind( umapCurveTab, pts )
}

fwrite( umapCurveTab, file = paste0( tablePath, "Trajectories_UMAP3D_Curves.txt" ), sep = "\t" )
