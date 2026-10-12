import csv
import json

# change to your project directory
home_directory = "path/to/project_directory/"

# data folder
dataPath = home_directory + "Data/HOG_RNAseq/"

# folder of the HOG data preparation (the fastp reports of Trimming_fastp.sh are in the subfolder fastp_reports)
outputPath = home_directory + "DataPreparation/HOG_RNAseq/"
tablePath = outputPath + "tables/"


# read overview of the HOG samples
# columns: sample, condition, replicate
with open(dataPath + "HOG_SampleOverview.txt") as f:
    samples = [row["sample"] for row in csv.DictReader(f, delimiter="\t")]

# number of reads before and after trimming (fastp report of each sample)
fastp_rows = []
for sample in samples:
    with open(outputPath + "fastp_reports/" + sample + "_fastp.json") as f:
        report = json.load(f)
    before = report["summary"]["before_filtering"]
    after = report["summary"]["after_filtering"]
    fastp_rows.append({"sample": sample, "fastp_version": report["summary"]["fastp_version"],
                       "reads_before": before["total_reads"], "bases_before": before["total_bases"], "q30_rate_before": before["q30_rate"],
                       "reads_after": after["total_reads"], "bases_after": after["total_bases"], "q30_rate_after": after["q30_rate"],
                       "duplication_rate": report["duplication"]["rate"], "insert_size_peak": report["insert_size"]["peak"]})

# number of read pairs that were mapped by Salmon (information that Salmon writes to each quantification)
salmon_rows = []
for sample in samples:
    with open(dataPath + "HOG_SalmonQuant/" + sample + "/aux_info/meta_info.json") as f:
        info = json.load(f)
    salmon_rows.append({"sample": sample, "salmon_version": info["salmon_version"], "library_type": ";".join(info["library_types"]),
                        "processed_pairs": info["num_processed"], "mapped_pairs": info["num_mapped"], "percent_mapped": info["percent_mapped"],
                        "index_seq_hash": info["index_seq_hash"], "index_name_hash": info["index_name_hash"]})

for name, rows in [("HOG_fastp_Overview.txt", fastp_rows), ("HOG_Salmon_Overview.txt", salmon_rows)]:
    with open(tablePath + name, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
