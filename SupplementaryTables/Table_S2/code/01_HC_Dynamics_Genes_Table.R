library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# results of the genes that change along pseudotime in HC (Figure 4)
dynPath = paste0( home_directory, "Figure_4/tables/HC_Dynamics/" )

# output folder
outputPath = paste0( home_directory, "SupplementaryTables/Table_S2/" )
tablePath = paste0( outputPath, "tables/" )


# Supplementary Table S2: genes that change along pseudotime in HC (association test, FDR < 0.05) with their gene module in Figure 4

# results of the association test (made by 08_HC_Dynamics_AssociationTest_Analysis.R of Figure 4)
geneTab = as.data.frame( fread( paste0( dynPath, "HC_Dynamics_AssociationTest_HC.txt.gz" ) ) )
# gene modules (made by 09_HC_Dynamics_Modules_Analysis.R of Figure 4)
# columns: raw_cluster, unit_id, trajectory_id, endpoint, condition, gene, ..., block_order, block, ...
moduleTab = as.data.frame( fread( paste0( dynPath, "HC_Dynamics_Modules_Genes.txt.gz" ) ) )

tab = geneTab[ geneTab$significant_TDE, ]
cand = match( paste( tab$unit_id, tab$gene ), paste( moduleTab$unit_id, moduleTab$gene ) )
# every significant gene has a module
stopifnot( !anyNA( cand ), nrow( moduleTab ) == nrow( tab ) )
tab$module = paste0( "C", moduleTab$block_order[ cand ] )

out = data.frame( lineage = tab$trajectory_id, endpoint = tab$endpoint, gene = tab$gene, pvalue = tab$pvalue, fdr = tab$padj_BH_lineage, module = tab$module )
out = out[ order( out$lineage, out$fdr, out$pvalue, out$gene, method = "radix" ), ]
table( out$lineage, out$module )

dir.create( tablePath, showWarnings = FALSE )
fwrite( out, file = paste0( tablePath, "Supplementary_Table_S2_HC_dynamic_genes.txt.gz" ), sep = "\t", na = "NA" )
