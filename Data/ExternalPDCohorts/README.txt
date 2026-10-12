Three published snRNA-seq data sets of PD and control donors (not part of this repository)

Data sets (downloaded from GEO, see Data/RawDataInfo.txt)
- Martirosyan et al., GEO series GSE243639
- Smajic et al., GEO series GSE157783
- Lee et al., GEO series GSE148434

The files of this folder are made from the downloaded data by the scripts in DataPreparation/ExternalPDCohorts/code, in this order:
1. Lee_DonorTable.R  ->  Lee_Donors.txt, the 12 donors of the Lee data set (gsm_id, sample_id, diagnosis, sex, age; needs Table S1 of the article in the folder of the molecule files, see Data/RawDataInfo.txt) (a few seconds)
2. Lee_Counts_QC_Scrublet.py  ->  Lee_Counts/ (count matrix of each donor), needs Python 3.10 with the packages of PythonPackages_Scrublet.txt (about 1 minute)
3. Lee_Clustering_AllCells.R  ->  cluster of each cell (about 5 minutes)
4. Lee_scCDC_Correction.R  ->  OPC_OL_Counts/Lee.rds (about 5 minutes)
5. Smajic_scCDC_Correction.R  ->  OPC_OL_Counts/Smajic.rds (about 5 minutes)
6. Martirosyan_scCDC_Correction.R  ->  OPC_OL_Counts/Martirosyan.rds (about 25 minutes)
7. Clustering_Pseudobulk_OPC_Oligodendrocytes.R  ->  Pseudobulk/<cohort>/<population>.rds, pseudobulk counts per donor (about 2 minutes)
