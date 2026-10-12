#!/bin/bash

# change to your project directory
home_directory="path/to/project_directory"

# trimmed reads and Salmon index
inputPath="${home_directory}/DataPreparation/HOG_RNAseq"

# output folder
outputPath="${home_directory}/Data/HOG_RNAseq/HOG_SalmonQuant"

mkdir -p "${outputPath}"

# transcript quantification of the six HOG samples with Salmon 1.10.3
# library type ISR: paired-end reads, stranded (read 1 from the reverse strand)
# a new run overwrites the quantification used for the article; with several threads Salmon gives slightly different numbers in each run
for sample in Ctrl1 Ctrl2 Ctrl3 asyn1 asyn2 asyn3
do
    echo "${sample}"

    salmon quant -i "${inputPath}/salmon_index" -l ISR \
           -1 "${inputPath}/TrimmedReads/${sample}_clean_R1.fastq.gz" -2 "${inputPath}/TrimmedReads/${sample}_clean_R2.fastq.gz" \
           --gcBias --seqBias --numBootstraps 30 --validateMappings -p 14 -o "${outputPath}/${sample}"
done
