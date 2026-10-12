# Transcriptional response of human oligodendroglial cells to extracellular alpha-synuclein, compared with oligodendroglia of MSA and PD patients


* This repository belongs to the publication
   * ***Extracellular alpha-synuclein induces transcriptional changes in human oligodendrocyte culture***
   * The journal reference will be added when the article is published.
* It contains all scripts and links to the data needed for the analyses, and the data shown in the figures. It shows step by step how to run the analyses with the scripts.
* It analyses RNA-seq of HOG cells (a human oligodendroglial cell line) exposed to alpha-synuclein pre-formed fibrils (PFF), GEO series [GSE349624](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349624).
* It re-analyses published single-nucleus RNA-seq data
   * striatal snRNA-seq of 10 healthy controls (HC), 7 MSA and 12 PD donors: Rydbirk et al., "Single-nucleus brain transcriptomics reveals microglia dysfunction in multiple system atrophy", __Nature Communications__ 17, 5234 (2026), https://doi.org/10.1038/s41467-026-71525-6
   * three snRNA-seq data sets of PD and control donors: Martirosyan et al. (GEO GSE243639), Smajić et al. (GEO GSE157783) and Lee et al. (GEO GSE148434)


## Repository structure

```
DataPreparation/          scripts from the raw data to the inputs of the figures
Data/                     inputs of the analyses (Data/DataFileInfo.txt), download instructions for the files that are not included
Figure_1/ ... Figure_5/   one folder per figure:
    code/                 scripts, NN_<Topic>_Analysis.R computes the results, NN_<Topic>_Plots.R draws the figure panels
    tables/, objects/     results
    figures/              figure panels (pdf); panels with many points have the points as image, *_FullVector.pdf is the same panel with all points as vector graphic
SupplementaryFigures/     Supplementary Figure 1 (SupplementaryFigure_1A, _1BC, _1DE), organized like the figure folders
SupplementaryTables/      Supplementary Tables 1 to 4 (Supplementary_Table_1.xlsx, ..., Supplementary_Table_4.xlsx) with the scripts and full results
Figure_Data/              the data shown in each figure panel, one file per panel
```


## Analysis pipeline

* The analyses are R scripts, a few steps are Python scripts (trajectory plots of Figure 4, reading of h5ad and Cell Ranger files) and shell scripts (fastp, Salmon).
* Each script reads the results of earlier scripts. Run the folders in this order and, within a folder, the scripts in the order of their numbers:
  1. DataPreparation (reference files, HOG RNA-seq, snRNA-seq of HC, MSA and PD donors, three published PD data sets)
  2. Figure_1 to Figure_5
  3. SupplementaryFigures, then SupplementaryTables
* Scripts without a number (saveFigure.R, panelLayout.R, plotTFActivity.R, readAuthorCounts.R, Lee_readCounts.R) are read by the other scripts and are not run on their own.

### The analysis has 2 entry points
1. Starting from the raw data
   * Trimming and quantification of the HOG reads, correction of ambient RNA in the snRNA-seq data, clustering and trajectory models
   * The data preparation takes several hours, fitting the tradeSeq models of Figure 4 and 5 takes several days
   * Go to "Prerequisites", then "Data preparation"
2. Starting from the results in this repository
   * The tables and objects in the figure folders are the results used for the article. The plotting scripts of Figure 1, 2, 3C, 3D, 3F, 3G, 4B–D and 5A–C, F, G draw the panels from them directly (Figure 1D–H and 5B, C, F, G also need the MSigDB gene sets of DataPreparation/ReferenceFiles/code/MSigDB_GeneSets.R)
   * The panels that show single nuclei or single donors (Figure 3A, 3B, 3E, 4A, 5D, 5E and Supplementary Figure 1) need the data preparation, because the single-nucleus data of the published studies are not part of this repository
   * Figure_Data holds the values of every panel as tables


## Prerequisites

### Software

* R 4.5.3 with Bioconductor 3.22, packages and versions in RPackages.txt
* Python 3.11 for the Python scripts (packages in PythonPackages.txt), Python 3.10 with Scrublet for DataPreparation/ExternalPDCohorts/code/Lee_Counts_QC_Scrublet.py (packages in PythonPackages_Scrublet.txt)
* fastp 1.3.2 and Salmon 1.10.3 for the HOG RNA-seq reads
* STREAM 1.1 for Figure 4A: download the source code into Data/Software (see Data/Software/README.txt)
* The R scripts run in parallel with forked processes and need macOS or Linux

### Data

The raw data are read from one folder (rawDataPath at the top of the scripts), Data/RawDataInfo.txt lists all files and where to get them:
* HOG_RNAseq: the paired-end reads of the six HOG samples, GEO [GSE349624](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349624)
* GENCODE_v47: transcript sequences, genome (primary assembly) and gene annotation of GENCODE release 47
* snRNAseq_HC_MSA_PD: adata.h5ad of the Single Cell Portal (SCP3252) and the donor table metadata.tsv of the data archive of Rydbirk et al. (see Data/snRNAseq_HC_MSA_PD/README.txt)
* ExternalPDCohorts: the files of the GEO series GSE243639, GSE157783 and GSE148434, and Table S1 of Lee et al. (age and sex of the donors)

Reference files that are downloaded into the Data folder (see the DOWNLOAD files there):
* Data/ReferenceFiles/GeneSets/Human.MitoPathways3.0.gmx: MitoCarta 3.0 pathways (Figure 2)
* Data/ReferenceFiles/TF_Networks/Oli_GRN.txt: gene regulatory network of oligodendrocytes of brainSCOPE (Figure 2, 3 and Supplementary Table 1)
* The MSigDB 2026.1 gene sets and the KEGG pathways are downloaded by the scripts of DataPreparation/ReferenceFiles (the KEGG API is for academic use only)

### Paths

* Set the project folder at the top of each script: home_directory, and where needed rawDataPath (raw data) and pythonPath (Python with anndata)


## Data preparation

### Reference files (DataPreparation/ReferenceFiles/code)

* ***GENCODE_v47_GeneTables.R***
   * transcript to gene table, protein-coding genes and the 13 protein-coding genes of the mitochondrial genome from the GENCODE annotation
   * ***saves*** Data/ReferenceFiles/GENCODE_v47/
* ***MSigDB_GeneSets.R***
   * GO (biological process, cellular component, molecular function), Reactome and hallmark gene sets of MSigDB 2026.1 with the R package msigdbr 26.1.0
   * ***saves*** Data/ReferenceFiles/GeneSets/MSigDB_2026.1_*
* ***KEGG_Pathways.R***
   * downloads the human KEGG pathways with clusterProfiler; KEGG has no version number, the article uses the download of 2026-09-05
   * ***saves*** Data/ReferenceFiles/GeneSets/KEGG_hsa_Pathways.rds

### RNA-seq of the HOG cells (DataPreparation/HOG_RNAseq/code)

* ***Trimming_fastp.sh***
   * adapter and quality trimming with fastp (quality 20, sliding window of 4 bases from both ends, poly-G trimming, minimum length 50)
* ***SalmonIndex_GENCODE_v47.sh***
   * decoy-aware Salmon index of the GENCODE v47 transcripts with the genome as decoy (k = 31)
* ***Quantification_Salmon.sh***
   * transcript quantification with Salmon (library type ISR, GC and sequence bias correction, 30 bootstraps)
   * ***saves*** Data/HOG_RNAseq/HOG_SalmonQuant/ (this overwrites the quantification of the article; Salmon with several threads gives slightly different numbers in each run)
* ***ReadNumbers_fastp_Salmon.py***
   * numbers of reads before and after trimming and mapping rates
   * ***saves*** DataPreparation/HOG_RNAseq/tables/

### snRNA-seq of HC, MSA and PD donors (DataPreparation/snRNAseq_HC_MSA_PD/code)

* ***CellAnnotation_AuthorLabels_UMAP.py***
   * cell type labels, donor, diagnosis and UMAP coordinates of the authors from adata.h5ad
   * ***saves*** Data/snRNAseq_HC_MSA_PD/NatComm_CellAnnotation.txt.gz
* ***DonorCovariates_AuthorMetadata.R***
   * donor table of the authors (diagnosis, sex, age, post-mortem interval, brain bank, batches)
   * ***saves*** Data/snRNAseq_HC_MSA_PD/NatComm_DonorCovariates.txt
* ***scCDC_ContaminationDetection.R***, ***scCDC_Correction.R***
   * detection of the genes contaminated by ambient RNA and correction of their counts with scCDC, separately for each donor (about 10 minutes each)
   * ***saves*** Data/snRNAseq_HC_MSA_PD/NatComm_CorrectedCounts_scCDC.rds

### Three published PD data sets (DataPreparation/ExternalPDCohorts/code)

* The order of the scripts and their run times are in Data/ExternalPDCohorts/README.txt
* ***Lee_DonorTable.R***, ***Lee_Counts_QC_Scrublet.py***, ***Lee_Clustering_AllCells.R***
   * count matrices of the Lee data set from the Cell Ranger molecule files, quality control as described by the authors (200 to 10,000 genes, at most 10% mitochondrial counts), doublet removal with Scrublet, clustering of all cells and cell types by marker genes
* ***Lee_scCDC_Correction.R***, ***Smajic_scCDC_Correction.R***, ***Martirosyan_scCDC_Correction.R***
   * correction of ambient RNA with scCDC, separately for each donor, then the OPCs and oligodendrocytes are kept
   * ***saves*** Data/ExternalPDCohorts/OPC_OL_Counts/
* ***Clustering_Pseudobulk_OPC_Oligodendrocytes.R***
   * clustering of the OPCs and oligodendrocytes of each data set with the same steps as for our own data (Figure 3), pseudobulk counts per donor and cell population
   * ***saves*** Data/ExternalPDCohorts/Pseudobulk/ (used by Supplementary Table 1)


## Figure 1: HOG cells exposed to alpha-synuclein PFF

* ***01_DifferentialExpression_Analysis.R***
   * gene counts from the Salmon quantification with tximport, genes with at least 10 reads in total
   * differential expression alpha-synuclein PFF vs control with DESeq2 (FDR < 0.05), PCA and sample distances for quality control
   * ***saves*** tables/DifferentialExpression/, objects/ (data of the volcano plot)
* ***02_DifferentialExpression_Plots.R***
   * ***produces*** Figure 1C (volcano plot of the protein-coding genes) and the quality control plots in figures/QC/
* ***03_Enrichment_Analysis.R***
   * over-representation of the up-regulated, down-regulated and all differentially expressed protein-coding genes with clusterProfiler, all tested genes as background
   * KEGG pathways, Reactome and GO (biological process, cellular component, molecular function), gene sets with at least 5 genes
   * ***saves*** tables/Enrichment/
* ***04_Enrichment_Plots.R***
   * ***produces*** Figure 1D–H (up to 8 gene sets with FDR < 0.05 per collection)


## Figure 2: mitochondrial and nuclear OXPHOS genes, transcription factors

* ***01_MitochondrialGenes_Analysis.R***, ***02_MitochondrialGenes_Plots.R***
   * z-scores of the 13 protein-coding genes of the mitochondrial genome and their log2 fold changes
   * ***produces*** Figure 2A
* ***03_GSEA_MitoCarta_Analysis.R***, ***04_GSEA_MitoCarta_Plots.R***
   * GSEA (fgsea) of the MitoCarta 3.0 pathways with all protein-coding genes ranked by the DESeq2 Wald statistic, once with and once without the genes of the mitochondrial genome
   * the leading edge genes of the OXPHOS subunits (nuclear genes only) are the signature of Figure 3E, the ranking of the nuclear genes is used in Figure 5
   * ***produces*** Figure 2B and 2C
* ***05_TF_Activity_Analysis.R***, ***06_TF_Activity_Plots.R***
   * the oligodendrocyte network of brainSCOPE is turned into a signed network (activating or repressing links)
   * transcription factor activities with VIPER (decoupleR) on the Wald statistics, at least 5 target genes, correction for shared targets
   * the signed network is read again by Figure 3 and Supplementary Table 1
   * ***produces*** Figure 2F


## Figure 3: oligodendroglia of HC, MSA and PD donors

* ***01_AllCells_Analysis.R***, ***02_AllCells_Plots.R***
   * numbers of cells and UMAPs of all cell types with the annotation of the authors (not shown in the article)
* ***03_OPC_OL_Clustering_Analysis.R***
   * OPCs and oligodendrocytes of the authors' annotation; a cell is removed if it highly expresses both marker genes of another cell type
   * 3000 variable genes without technical genes (Data/ReferenceFiles/GeneLists/TechnicalGenes.txt), integration of the donors with harmony, clustering with 30 components and resolution 0.3
   * clusters with at least 200 cells and at least 3 donors with 100 or more cells are kept: one OPC cluster and six oligodendrocyte clusters (OL1 = largest)
   * ***saves*** objects/ and tables/OPC_OL_Clustering/
* ***04_OPC_OL_UMAP_Plots.R***
   * ***produces*** Figure 3A and 3B
* ***05_CellStateProportions_Analysis.R***, ***06_CellStateProportions_Plots.R***
   * proportion of each cell state per donor, test of each disease against HC with propeller adjusted for age
   * ***produces*** Figure 3C
* ***07_DonorCovariates_Balance_Analysis.R***
   * balance of age, post-mortem interval, sex, brain bank and batches between the disease groups and HC for each cell state
* ***08_Pseudobulk_DESeq2_Analysis.R***, ***09_Pseudobulk_DESeq2_Plots.R***
   * pseudobulk counts per donor and cell state (donors with at least 10 cells), DESeq2 with age and diagnosis, genes detected in at least half of the donors
   * differentially expressed: FDR < 0.05 and absolute log2 fold change above 0.3
   * ***produces*** Figure 3D
* ***10_AUCell_OXPHOS_Analysis.R***, ***11_AUCell_OXPHOS_Plots.R***
   * AUCell score of the nuclear OXPHOS signature of the HOG cells (Figure 2C) in each nucleus
   * ***produces*** Figure 3E
* ***12_GSEA_Hallmark_Analysis.R***, ***13_GSEA_Hallmark_Plots.R***
   * GSEA of the hallmark gene sets with the protein-coding genes ranked by the DESeq2 Wald statistic, for each cell state and comparison
   * ***produces*** Figure 3F (the 4 most significant gene sets of each cell state and comparison)
* ***14_TF_Activity_OPC_OL_Analysis.R***, ***15_TF_Activity_OPC_OL_Plots.R***
   * transcription factor activities as for the HOG cells (Figure 2F), for each cell state, MSA vs HC and PD vs HC
   * ***produces*** Figure 3G (the 10 transcription factors that are most consistently significant over the cell states in PD)


## Figure 4: trajectories of the oligodendroglia and expression along pseudotime in HC

* ***01_Trajectories_Slingshot_Analysis.R***
   * trajectories with slingshot in the first three harmony components, starting in the OPC cluster: three lineages ending in OL5, OL6 and OL3
   * three-dimensional UMAP with the curves of the lineages
* ***02_Trajectories_STREAM_Analysis.py***, ***03_Trajectories_STREAM_Plots.py***, ***04_Trajectories_UMAP3D_Plots.py***
   * the slingshot lineages drawn as a tree with the plot function of STREAM, and in the three-dimensional UMAP
   * ***produces*** Figure 4A
* ***05_evaluateK_Analysis.R***, ***06_evaluateK_Plots.R***
   * number of knots of the tradeSeq models from evaluateK (3 to 10 knots, 100 genes, each lineage); the AIC values are read from Data/SavedModels/evaluateK, the call that made them is kept as a comment
   * 6 knots are used
* ***07_fitGAM_Analysis.R***
   * tradeSeq models (fitGAM) of all genes detected in at least 100 cells, three lineages with the conditions HC, MSA and PD, offset from TMM-normalized library sizes
   * runs several days; the models of the article were fitted on x86_64, on arm64 (for example Apple silicon) the same seed gives other fits
   * ***saves*** objects/fitGAM_HC_MSA_PD_ThreeLineages/
* ***08_HC_Dynamics_AssociationTest_Analysis.R***
   * association test of tradeSeq for HC (log2 fold change of 1 along each lineage, FDR < 0.05 within each lineage)
   * fitted curves of the significant genes in 50 bins of each lineage with the same number of HC cells, taken at the median pseudotime of each bin, z-score of each gene over the bins
* ***09_HC_Dynamics_Modules_Analysis.R***, ***10_HC_Dynamics_Heatmaps_Plots.R***
   * four modules per lineage by hierarchical clustering (Ward) of the z-scored curves, GO biological process enrichment of the modules (MSigDB 2026.1 GO BP gene sets), distribution of the cell states over the bins
   * ***produces*** Figure 4B–D


## Figure 5: MSA and PD compared with HC along pseudotime

* ***01_TwoLineages_fitGAM_Analysis.R***
   * tradeSeq models of lineage 1 (to OL5) and lineage 3 (to OL3) for MSA + HC and for PD + HC, 6 knots; runs several days
   * the comparison is restricted to the pseudotime range of the oligodendrocytes
   * ***saves*** objects/fitGAM_MSA_vs_HC_TwoLineages/, objects/fitGAM_PD_vs_HC_TwoLineages/
* ***02_Disease_Dynamics_ConditionTest_Analysis.R***
   * condition test of tradeSeq (log2 fold change of 1, FDR < 0.05 within each disease and lineage)
   * fitted HC and disease curves of the significant genes in 50 groups with the same number of cells (HC and disease together) within the compared pseudotime range, taken at the median pseudotime of each group; z-score of the log expression of each gene over both curves
* ***03_Disease_Dynamics_Modules_Analysis.R***, ***04_Disease_Dynamics_Modules_Plots.R***
   * modules by hierarchical clustering (Ward) of the z-scored HC and disease curves and dynamic tree cut (at least 20 genes)
   * GSEA of each module in the ranking of the nuclear genes of the HOG cells (Figure 2), GO biological process enrichment of each module
   * ***produces*** Figure 5A, 5B, 5C, 5F and 5G (figures/Figure_5_PanelFiles.txt lists which file is which panel)
* ***05_SelectedGenes_OXPHOS_PD_Analysis.R***, ***06_SelectedGenes_OXPHOS_PD_Plots.R***
   * fitted expression of COX6B1 and ATP5F1E for PD vs HC, with half of the nuclei shown as points
   * ***produces*** Figure 5D and 5E


## Supplementary Figure 1

* ***SupplementaryFigure_1A/code/01_CellStateProportions_PerDonor_Plots.R***
   * proportions of the cell states of every donor with the adjusted p-values of the propeller tests (Figure 3)
   * ***produces*** Supplementary Figure 1A
* ***SupplementaryFigure_1BC/code/***, ***SupplementaryFigure_1DE/code/***
   * fitted expression for MSA vs HC of COX6B1 and ATP5F1E (1B, 1C) and of the protein-folding genes HSPA1A and HSPA4L (1D, 1E)
   * ***produces*** Supplementary Figure 1B–E


## Supplementary Tables

* ***Table_S1/code/***
   * transcription factor activities PD vs HC in the three published data sets, with the same steps as for our own cell states (pseudobulk per donor, DESeq2 with age, VIPER with the brainSCOPE network)
   * ***produces*** Table_S1/Supplementary_Table_1.xlsx (selected transcription factors with FDR < 0.05) and the full results in Table_S1/tables/
* ***Table_S2/code/***
   * genes that change along pseudotime in HC (association test of Figure 4) with their gene module
   * ***produces*** Table_S2/Supplementary_Table_2.xlsx (all genes with FDR < 0.05) and the full results in Table_S2/tables/
* ***Table_S3/code/***
   * genes that differ between disease and HC along the oligodendrocyte part of lineages 1 and 3 (condition test of Figure 5) with their gene module
   * ***produces*** Table_S3/Supplementary_Table_3.xlsx (all genes with FDR < 0.05) and the full results in Table_S3/tables/
* ***Table_S4/code/***
   * GO biological process enrichment of the modules of Figure 5
   * ***produces*** Table_S4/Supplementary_Table_4.xlsx (selected processes) and the full results in Table_S4/tables/


## Figure data

* Figure_Data holds the data of each panel as tab-separated tables, only what is drawn in the panel, see Figure_Data/README.txt


## Citation

The citation will be added when the article is published.
