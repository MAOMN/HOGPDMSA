import os

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scrublet as scr

# change to your project directory
home_directory = "path/to/project_directory/"

# folder with the downloaded raw data (not part of this repository, see Data/RawDataInfo.txt)
rawDataPath = "path/to/raw_data/"

# data folder
dataPath = home_directory + "Data/"

# output folder
outputPath = home_directory + "DataPreparation/ExternalPDCohorts/"
tablePath = outputPath + "tables/"


# snRNA-seq data set of Lee et al. (substantia nigra, PD and control donors), GEO series GSE148434
# for each donor only the molecule file of Cell Ranger is available (molecule_info.h5), the count matrices are made from these files
geoPath = rawDataPath + "ExternalPDCohorts/Lee_GSE148434/"

# read donor table
# columns: gsm_id, sample_id, diagnosis, sex, age
donors = pd.read_csv(dataPath + "ExternalPDCohorts/Lee_Donors.txt", sep="\t")
print(donors)

os.makedirs(dataPath + "ExternalPDCohorts/Lee_Counts", exist_ok=True)

qc_tables = []
for gsm, sample in zip(donors["gsm_id"], donors["sample_id"]):
    with h5py.File(geoPath + gsm + "_" + sample + "_snRNA.molecule_info.h5", "r") as f:
        gene_ids = f["features/id"][()].astype(str)
        gene_names = f["features/name"][()].astype(str)

        # barcodes that were called as cells by Cell Ranger (first column: index of the barcode)
        cell_index = f["barcode_info/pass_filter"][:, 0].astype(np.int64)
        barcodes = f["barcodes"][cell_index].astype(str)

        # one row per molecule (UMI), each molecule of a cell is counted once for its gene
        molecule_barcode = f["barcode_idx"][()].astype(np.int64)
        molecule_gene = f["feature_idx"][()]
        position = np.searchsorted(cell_index, molecule_barcode)
        position[position == len(cell_index)] = 0
        in_cell = (cell_index[position] == molecule_barcode) & (f["count"][()] > 0)

    counts = sp.coo_matrix((np.ones(in_cell.sum(), dtype=np.int32), (molecule_gene[in_cell], position[in_cell])),
                           shape=(len(gene_names), len(cell_index))).tocsc()
    counts.sum_duplicates()
    counts.sort_indices()

    # quality control as described by the authors: 200 to 10000 detected genes, at most 10% mitochondrial counts
    total_umi = np.asarray(counts.sum(axis=0)).ravel()
    detected_genes = np.asarray((counts > 0).sum(axis=0)).ravel()
    is_mt = np.char.startswith(gene_names, "MT-")
    percent_mt = 100 * np.asarray(counts[is_mt, :].sum(axis=0)).ravel() / total_umi
    qc_pass = (detected_genes >= 200) & (detected_genes <= 10000) & (percent_mt <= 10)

    # doublet scores with Scrublet (all cells of the donor), cells with a score of 0.4 or more are removed
    scrub = scr.Scrublet(counts.T.tocsc(), expected_doublet_rate=0.1, random_state=0)
    doublet_score, _ = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30, verbose=False)
    keep = qc_pass & (doublet_score < 0.4)
    print(sample, "cells of Cell Ranger:", len(keep), "kept:", keep.sum())

    cell_ids = np.array([sample + "__" + b for b in barcodes])
    qc_tables.append(pd.DataFrame({"cell_id": cell_ids, "sample_id": sample, "total_umi": total_umi, "detected_genes": detected_genes,
                                   "percent_mt": percent_mt, "doublet_score": doublet_score, "kept": keep}))

    # counts of the kept cells (genes x cells, compressed sparse column format)
    kept_counts = counts[:, keep].tocsc()
    with h5py.File(dataPath + "ExternalPDCohorts/Lee_Counts/" + sample + "_Counts.h5", "w") as out:
        out.create_dataset("matrix/data", data=kept_counts.data, compression="gzip")
        out.create_dataset("matrix/indices", data=kept_counts.indices, compression="gzip")
        out.create_dataset("matrix/indptr", data=kept_counts.indptr, compression="gzip")
        out.create_dataset("matrix/shape", data=np.array(kept_counts.shape))
        out.create_dataset("features/id", data=gene_ids.astype("S"))
        out.create_dataset("features/name", data=gene_names.astype("S"))
        out.create_dataset("cells/cell_id", data=cell_ids[keep].astype("S"))

# quality measures of all cells
qc = pd.concat(qc_tables)
print(qc.groupby("sample_id")["kept"].agg(["size", "sum"]))
qc.to_csv(tablePath + "Lee_CellQC_Scrublet.txt.gz", sep="\t", index=False, compression="gzip")
