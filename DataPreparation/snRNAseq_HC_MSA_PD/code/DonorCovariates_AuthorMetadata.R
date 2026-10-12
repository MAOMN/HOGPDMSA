library(data.table)

# change to your project directory
home_directory <- "path/to/project_directory/"

# folder with the downloaded raw data (not part of this repository, see Data/snRNAseq_HC_MSA_PD/README.txt)
rawDataPath <- "path/to/raw_data/"

# data folder
dataPath = paste0( home_directory, "Data/snRNAseq_HC_MSA_PD/" )


# read the donor table of the authors (metadata.tsv of their data archive)
# columns used: sample, condition, sex, age, pmi, brain_bank, extract, flowcell, disease_duration, subtype
authorTab = fread( paste0( rawDataPath, "snRNAseq_HC_MSA_PD/metadata.tsv" ) )
table( authorTab$condition )

# the controls are called CTRL by the authors, here HC
donor_info = data.frame( sample_id = as.character( authorTab$sample ), Diagnosis_author = as.character( authorTab$condition ),
                         Diagnosis = ifelse( authorTab$condition == "CTRL", "HC", as.character( authorTab$condition ) ),
                         Sex = as.character( authorTab$sex ), Age = as.numeric( authorTab$age ), PMI = as.numeric( authorTab$pmi ),
                         brain_bank = as.character( authorTab$brain_bank ), extract = as.character( authorTab$extract ), flowcell = as.character( authorTab$flowcell ),
                         disease_duration = as.numeric( authorTab$disease_duration ), subtype = as.character( authorTab$subtype ) )

# one row per donor
donor_info = donor_info[ !duplicated( donor_info$sample_id ), ]
nrow( donor_info )

fwrite( donor_info, file = paste0( dataPath, "NatComm_DonorCovariates.txt" ), sep = "\t" )
