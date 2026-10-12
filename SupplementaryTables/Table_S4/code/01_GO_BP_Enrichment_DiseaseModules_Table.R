library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# results of the comparison of disease and HC along the trajectories (Figure 5)
dynPath = paste0( home_directory, "Figure_5/tables/Modules/" )

# output folder
outputPath = paste0( home_directory, "SupplementaryTables/Table_S4/" )
tablePath = paste0( outputPath, "tables/" )


# Supplementary Table S4: GO biological process enrichment of all gene modules of Figure 5 (all tested terms, not only the significant ones)

# GO enrichment results of the modules (made by 03_Disease_Dynamics_Modules_Analysis.R of Figure 5)
goTab = as.data.frame( fread( paste0( dynPath, "Disease_Dynamics_Modules_GO_BP.txt.gz" ) ) )

# overview of the modules
# columns: unit_id, module_id, module, genes, ..., unit_order, module_order, display_cluster, panel_order, display_unit
moduleTab = as.data.frame( fread( paste0( dynPath, "Disease_Dynamics_Modules_Overview.txt.gz" ) ) )

cand = match( goTab$module_id, moduleTab$module_id )
goTab$module_label = moduleTab$display_cluster[ cand ]
goTab$module_genes_total = moduleTab$genes[ cand ]
goTab$unit_order = moduleTab$unit_order[ cand ]
goTab$module_number = as.integer( sub( "^C", "", goTab$module_label ) )

goTab$contrast = paste0( goTab$disease, "_vs_HC" )
# GOBP_PROTEIN_FOLDING -> protein folding
goTab$gene_set_label = gsub( "_", " ", tolower( sub( "^GOBP_", "", goTab$term_name ) ) )
# size of the gene set within the gene universe: first number of the background ratio
goTab$gene_set_size_in_universe = as.integer( sapply( strsplit( goTab$BgRatio, "/", fixed = TRUE ), function( x ) x[ 1 ] ) )

goTab = goTab[ order( goTab$unit_order, goTab$module_number, goTab$padj_within_module_BH, goTab$pvalue, goTab$term_id, method = "radix" ), ]

tableS4 = data.frame( contrast = goTab$contrast, disease = goTab$disease, lineage = goTab$lineage_id, module_label = goTab$module_label, module_id = goTab$module_id,
                      msigdb_gene_set_id = goTab$term_id, msigdb_gene_set_name = goTab$term_name, gene_set_label = goTab$gene_set_label,
                      module_genes_total = goTab$module_genes_total, module_genes_in_universe = goTab$query_genes_mapped, background_genes_in_universe = goTab$background_mapped,
                      gene_set_size_in_universe = goTab$gene_set_size_in_universe, overlap_genes_n = goTab$overlap_n, overlap_genes = goTab$overlap_genes,
                      p_value = goTab$pvalue, fdr_bh_within_module = goTab$padj_within_module_BH, significant_fdr_lt_0_05 = goTab$FDR05 )

fwrite( tableS4, file = paste0( tablePath, "Supplementary_Table_S4_GO_BP_enrichment.txt.gz" ), sep = "\t", na = "NA", quote = FALSE )

# number of tested and significant terms per module
moduleSummary = unique( goTab[ , c( "contrast", "disease", "lineage_id", "module_label", "module_id", "module_genes_total", "query_genes_mapped", "background_mapped", "correction_family" ) ] )
colnames( moduleSummary ) = c( "contrast", "disease", "lineage", "module_label", "module_id", "module_genes_total", "module_genes_in_universe", "background_genes_in_universe", "correction_family" )
moduleSummary$tested_terms_with_overlap = sapply( moduleSummary$module_id, function( m ) sum( goTab$module_id == m ) )
moduleSummary$significant_terms_fdr_lt_0_05 = sapply( moduleSummary$module_id, function( m ) sum( goTab$FDR05[ goTab$module_id == m ] ) )

fwrite( moduleSummary, file = paste0( tablePath, "Supplementary_Table_S4_module_summary.txt" ), sep = "\t", na = "NA", quote = FALSE )
