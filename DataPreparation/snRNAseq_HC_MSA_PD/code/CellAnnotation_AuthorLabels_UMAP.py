import anndata as ad
import pandas as pd

# change to your project directory
home_directory = "path/to/project_directory/"

# folder with the downloaded raw data (not part of this repository, see Data/snRNAseq_HC_MSA_PD/README.txt)
rawDataPath = "path/to/raw_data/"

# data folder
dataPath = home_directory + "Data/snRNAseq_HC_MSA_PD/"


# read the snRNA-seq data set of the authors, only the cell annotation and the UMAP coordinates are needed here
adata = ad.read_h5ad(rawDataPath + "snRNAseq_HC_MSA_PD/adata.h5ad", backed="r")
obs = adata.obs
print(obs["disease__ontology_label"].value_counts())

# donor, diagnosis, cell type labels (major and minor) and UMAP of the authors for each cell
diagnosis_names = {"normal": "HC", "multiple system atrophy": "MSA", "Parkinson disease": "PD"}
cells = pd.DataFrame({"cell_id": list(adata.obs_names),
                      "sample_id": obs["donor_id"].astype(str).values,
                      "diagnosis": obs["disease__ontology_label"].astype(str).map(diagnosis_names).values,
                      "labels_major": obs["labels_major"].astype(str).values,
                      "labels_minor": obs["labels_minor"].astype(str).values,
                      "UMAP_1": adata.obsm["X_umap"][:, 0],
                      "UMAP_2": adata.obsm["X_umap"][:, 1]})
print(cells.shape)

cells.to_csv(dataPath + "NatComm_CellAnnotation.txt.gz", sep="\t", index=False, compression="gzip")
adata.file.close()
