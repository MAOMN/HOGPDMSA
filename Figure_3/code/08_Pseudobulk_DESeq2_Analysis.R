library(Seurat)
library(Matrix)
library(DESeq2)
library(org.Hs.eg.db)
library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

dataPath = paste0( home_directory, "Data/" )
outputPath = paste0( home_directory, "Figure_3/" )
tablePath = paste0( outputPath, "tables/Pseudobulk_DESeq2/" )
objectPath = paste0( outputPath, "objects/" )


# Seurat object with the OPCs and oligodendrocytes of the kept clusters (03_OPC_OL_Clustering_Analysis.R)
obj <- readRDS( paste0( objectPath, "OPC_OL_Seurat_KeptClusters.rds" ) )
counts <- GetAssayData( obj, layer = "counts" )
meta <- obj@meta.data

# donor information (columns: sample_id, Diagnosis_author, Diagnosis, Sex, Age, PMI, brain_bank, extract, flowcell, disease_duration, subtype)
# made by DataPreparation/snRNAseq_HC_MSA_PD/code/DonorCovariates_AuthorMetadata.R from the data of Rydbirk et al. (2026)
donor_info = as.data.frame( fread( paste0( dataPath, "snRNAseq_HC_MSA_PD/NatComm_DonorCovariates.txt" ) ) )
rownames( donor_info ) <- donor_info$sample_id

# technical genes (ribosomal and mitochondrial genes, see Data/DataFileInfo.txt; columns: gene_symbol, ensembl_id)
technicalGenes = as.data.frame( fread( paste0( dataPath, "ReferenceFiles/GeneLists/TechnicalGenes.txt" ) ) )

states = c( "OPC", paste0( "OL", 1:6 ) )
comparisons = c( MSA_vs_HC = "MSA", PD_vs_HC = "PD" )


#
### Pseudobulk: sum of the counts of all cells of a donor, separately for each cell state
#

pseudobulk = list()
cellsPerDonor = list()
overview = c()
for( state in states )
{
    pb = matrix( 0, nrow = nrow( counts ), ncol = nrow( donor_info ), dimnames = list( rownames( counts ), donor_info$sample_id ) )
    n = rep( 0, nrow( donor_info ) )
    names( n ) = donor_info$sample_id
    for( donor in donor_info$sample_id )
    {
        cand = which( meta$final_cluster_label == state & meta$sample_id == donor )
        n[ donor ] = length( cand )
        if( length( cand ) > 0 ) pb[ , donor ] = rowSums( counts[ , cand, drop = FALSE ] )
    }
    pseudobulk[[ state ]] = pb
    cellsPerDonor[[ state ]] = n
    overview = rbind( overview, data.frame( natural_label = state, n_cells = sum( n ), n_samples_with_cells = sum( n > 0 ), min_cells_per_sample = min( n ), max_cells_per_sample = max( n ) ) )
}
fwrite( overview, file = paste0( tablePath, "Pseudobulk_CellNumbers.txt" ), sep = "\t" )
rm( obj, counts )
gc()


#
### Differential expression per cell state: MSA vs HC and PD vs HC
#

# model: diagnosis and age of the donors (age standardized over the donors used for the cell state)
# donors with less than 10 cells of a cell state are not used

# gene type (protein-coding or not) from the gene symbol
geneType <- mapIds( org.Hs.eg.db, keys = rownames( pseudobulk[[ 1 ]] ), column = "GENETYPE", keytype = "SYMBOL", multiVals = "first" )

deseqResults = list()
summaryTab = c()
for( state in states )
{
    for( comp in names( comparisons ) )
    {
        disease = comparisons[ comp ]

        d = donor_info[ donor_info$Diagnosis %in% c( "HC", disease ) & cellsPerDonor[[ state ]][ donor_info$sample_id ] >= 10, ]
        d = d[ order( d$sample_id, method = "radix" ), ]
        d$Diagnosis = factor( d$Diagnosis, levels = c( "HC", disease ) )
        d$Age_z = ( d$Age - mean( d$Age ) ) / sd( d$Age )
        n_cells = cellsPerDonor[[ state ]][ d$sample_id ]

        pb = pseudobulk[[ state ]][ , d$sample_id ]
        # gene has to be detected in at least half of the donors
        keep = rowSums( pb > 0 ) >= ceiling( ncol( pb ) / 2 ) & rowSums( pb ) >= 10
        pb = pb[ keep, ]

        dds <- DESeqDataSetFromMatrix( countData = pb, colData = d[ , c( "Age_z", "Diagnosis" ) ], design = ~ Age_z + Diagnosis )
        dds <- DESeq( dds, fitType = "local", quiet = TRUE )
        res <- results( dds, contrast = c( "Diagnosis", disease, "HC" ), alpha = 0.05 )

        deseqResults[[ paste0( comp, "__", state ) ]] = res

        resTable = data.frame( natural_label = state, contrast = comp, case = disease, reference = "HC", model_formula = "~ Age_z + Diagnosis",
                               gene = rownames( res ), gene_type = geneType[ rownames( res ) ], as.data.frame( res ), row.names = NULL )
        resTable$is_protein_coding = resTable$gene_type == "protein-coding"

        # technical genes are not reported, also not genes without adjusted p-value
        isTechnical = toupper( resTable$gene ) %in% technicalGenes$gene_symbol
        resTable = resTable[ !isTechnical & !is.na( resTable$padj ), ]

        # differentially expressed: adjusted p-value below 0.05 and absolute log2 fold change above 0.3
        resTable$lfc_display_threshold = 0.3
        resTable$significance_lfc03 = ifelse( resTable$padj < 0.05 & resTable$log2FoldChange > 0.3, "Up", ifelse( resTable$padj < 0.05 & resTable$log2FoldChange < -0.3, "Down", "NS" ) )
        resTable = resTable[ , c( "natural_label", "contrast", "case", "reference", "model_formula", "gene", "gene_type", "is_protein_coding",
                                  "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "lfc_display_threshold", "significance_lfc03" ) ]

        fwrite( resTable, file = paste0( tablePath, "Pseudobulk_DESeq2_", comp, "_", state, ".txt.gz" ), sep = "\t" )

        sigPC = resTable[ which( resTable$is_protein_coding & resTable$significance_lfc03 != "NS" ), ]
        summaryTab = rbind( summaryTab, data.frame( natural_label = state, contrast = comp, model_formula = "~ Age_z + Diagnosis",
                                                    n_HC = sum( d$Diagnosis == "HC" ), n_case = sum( d$Diagnosis == disease ), n_samples = nrow( d ),
                                                    n_genes_after_count_filter = sum( keep ),
                                                    n_genes_after_snrna_filter_padj_nonNA = nrow( resTable ),
                                                    n_protein_coding_display = sum( resTable$is_protein_coding, na.rm = TRUE ),
                                                    n_sig_lfc03_protein_coding = nrow( sigPC ),
                                                    n_up_lfc03_protein_coding = sum( sigPC$log2FoldChange > 0 ),
                                                    n_down_lfc03_protein_coding = sum( sigPC$log2FoldChange < 0 ),
                                                    min_cells_per_sample = min( n_cells ), median_cells_per_sample = median( n_cells ) ) )
    }
}
fwrite( summaryTab, file = paste0( tablePath, "Pseudobulk_DESeq2_Overview.txt" ), sep = "\t" )

# DESeq2 results of all genes of each cell state and comparison as R object
saveRDS( deseqResults, file = paste0( objectPath, "Pseudobulk_DESeq2_Results.rds" ) )
