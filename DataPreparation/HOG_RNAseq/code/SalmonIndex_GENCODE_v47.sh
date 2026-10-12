#!/bin/bash

# change to your project directory
home_directory="path/to/project_directory"

# folder with the downloaded raw data (not part of this repository, see Data/RawDataInfo.txt)
rawDataPath="path/to/raw_data"

# output folder
outputPath="${home_directory}/DataPreparation/HOG_RNAseq"

mkdir -p "${outputPath}"

# GENCODE release 47: transcript sequences and genome (primary assembly)
transcripts="${rawDataPath}/GENCODE_v47/gencode.v47.transcripts.fa.gz"
genome="${rawDataPath}/GENCODE_v47/GRCh38.primary_assembly.genome.fa.gz"

# decoy-aware index: the genome sequences are added behind the transcripts as decoys
# names of the genome sequences (first word of the header lines)
gunzip -c "${genome}" | grep '^>' | sed 's/^>//' | cut -d' ' -f1 > "${outputPath}/decoys.txt"
cat "${transcripts}" "${genome}" > "${outputPath}/gentrome.fa.gz"

# Salmon 1.10.3, --gencode shortens the transcript names to the transcript id
salmon index -t "${outputPath}/gentrome.fa.gz" -d "${outputPath}/decoys.txt" -i "${outputPath}/salmon_index" -k 31 --gencode -p 14
