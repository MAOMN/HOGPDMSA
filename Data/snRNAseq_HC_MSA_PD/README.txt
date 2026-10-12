snRNA-seq data set of HC, MSA and PD donors

Source
- Rydbirk R, et al. Single-nucleus brain transcriptomics reveals microglia dysfunction in multiple system atrophy. Nature Communications 17, 5234 (2026). https://doi.org/10.1038/s41467-026-71525-6
- processed snRNA-seq data of the authors: Broad Institute Single Cell Portal, accession SCP3252 (https://singlecell.broadinstitute.org/single_cell/study/SCP3252)
- donor table of the authors: data archive that is linked in https://github.com/rrydbirk/MSAvsPD ("Remaining data", https://osf.io/bgrqa)
- raw sequencing data of the authors: European Genome-phenome Archive, accession EGAS50000001406 (restricted access, not used here)

Files of this folder
The data of the authors are not distributed with this repository. The three files below are made from the original data by the scripts in DataPreparation/snRNAseq_HC_MSA_PD/code
and are read by the scripts of Figure 3 (and the scripts that use its results).
- NatComm_CellAnnotation.txt.gz: cell identifier, donor, diagnosis, cell type labels and UMAP coordinates of the authors (from adata.h5ad of SCP3252)
- NatComm_DonorCovariates.txt: donor table (diagnosis, sex, age, post-mortem interval, brain bank, RNA extraction batch, flowcell) from metadata.tsv of the data archive of the authors
- NatComm_CorrectedCounts_scCDC.rds: raw counts of adata.h5ad (genes x cells) after correction of ambient RNA contamination with scCDC

How to make the files
1. Download adata.h5ad from SCP3252 and metadata.tsv from the data archive of the authors (see above) into the folder <rawDataPath>/snRNAseq_HC_MSA_PD/ (rawDataPath is set at the top of the scripts)
   (SCP3252 offers the file as "adata (1).h5ad", rename it to adata.h5ad)
2. Run the scripts of DataPreparation/snRNAseq_HC_MSA_PD/code
   - CellAnnotation_AuthorLabels_UMAP.py  ->  NatComm_CellAnnotation.txt.gz
   - DonorCovariates_AuthorMetadata.R  ->  NatComm_DonorCovariates.txt
   - scCDC_ContaminationDetection.R  ->  contaminating genes of each donor (about 10 minutes)
   - scCDC_Correction.R  ->  NatComm_CorrectedCounts_scCDC.rds (about 10 minutes), it reads the contaminating genes of the previous script

Please cite the article above when these files are used.
