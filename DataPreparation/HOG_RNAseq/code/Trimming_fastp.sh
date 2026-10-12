#!/bin/bash

# change to your project directory
home_directory="path/to/project_directory"

# folder with the downloaded raw data (not part of this repository, see Data/RawDataInfo.txt)
rawDataPath="path/to/raw_data"

# output folder
outputPath="${home_directory}/DataPreparation/HOG_RNAseq"

mkdir -p "${outputPath}/TrimmedReads" "${outputPath}/fastp_reports"

# trimming of the paired-end reads of the six HOG samples with fastp 1.3.2
# the raw files are called <sample>_S<number>_R1_001.fastq.gz and <sample>_S<number>_R2_001.fastq.gz
for s in Ctrl1_S1 Ctrl2_S2 Ctrl3_S3 asyn1_S4 asyn2_S5 asyn3_S6
do
    sample=${s%_S*}
    echo "${sample}"

    fastp -i "${rawDataPath}/HOG_RNAseq/${s}_R1_001.fastq.gz" -I "${rawDataPath}/HOG_RNAseq/${s}_R2_001.fastq.gz" \
          -o "${outputPath}/TrimmedReads/${sample}_clean_R1.fastq.gz" -O "${outputPath}/TrimmedReads/${sample}_clean_R2.fastq.gz" \
          --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50 --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --trim_poly_g \
          --thread 14 --json "${outputPath}/fastp_reports/${sample}_fastp.json" --html "${outputPath}/fastp_reports/${sample}_fastp.html"
done
