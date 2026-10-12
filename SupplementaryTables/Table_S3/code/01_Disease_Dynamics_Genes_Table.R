library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# results of the comparison of disease and HC along the trajectories (Figure 5)
dynPath = paste0( home_directory, "Figure_5/" )

# output folder
outputPath = paste0( home_directory, "SupplementaryTables/Table_S3/" )
tablePath = paste0( outputPath, "tables/" )


# Supplementary Table S3: genes that differ between disease and HC along the oligodendrocyte part of lineages 1 and 3 (condition test, FDR < 0.05), with their gene module

# results of the condition test (made by 02_Disease_Dynamics_ConditionTest_Analysis.R of Figure 5)
condTab <- readRDS( paste0( dynPath, "objects/Disease_Dynamics_ConditionTest.rds" ) )
# gene modules and their names in the figures (made by 03_Disease_Dynamics_Modules_Analysis.R of Figure 5)
moduleGenes = as.data.frame( fread( paste0( dynPath, "tables/Modules/Disease_Dynamics_Modules_Genes.txt.gz" ) ) )
moduleTab = as.data.frame( fread( paste0( dynPath, "tables/Modules/Disease_Dynamics_Modules_Overview.txt.gz" ) ) )

tab = condTab[ condTab$XDE_FDR05, ]
cand = match( paste( tab$unit_id, tab$gene ), paste( moduleGenes$unit_id, moduleGenes$gene ) )
stopifnot( !anyNA( cand ) )
moduleId = ifelse( moduleGenes$assigned[ cand ], moduleGenes$module_id[ cand ], NA )
tab$module = moduleTab$display_cluster[ match( moduleId, moduleTab$module_id ) ]

out = data.frame( comparison = paste( tab$disease, "vs HC" ), lineage = tab$lineage_id, endpoint = tab$endpoint, gene = tab$gene,
                  pvalue = tab$pvalue, fdr = tab$padj_BH, module = tab$module )
out = out[ order( out$comparison, out$lineage, out$fdr, out$pvalue, out$gene, method = "radix" ), ]
table( out$comparison, out$lineage )
table( is.na( out$module ) )

dir.create( tablePath, showWarnings = FALSE )
fwrite( out, file = paste0( tablePath, "Supplementary_Table_S3_disease_dynamic_genes.txt.gz" ), sep = "\t", na = "NA" )
