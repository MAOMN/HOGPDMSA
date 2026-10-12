import os
import pickle
import subprocess
import sys

import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd

# change to your project directory
home_directory = "path/to/project_directory/"

# data folder
dataPath = home_directory + "Data/"

# output folder
outputPath = home_directory + "Figure_4/"
tablePath = outputPath + "tables/Trajectories/"
objectPath = outputPath + "objects/"

# STREAM 1.1 is not installed as a package, its source code is loaded from the data folder (download see Data/Software/README.txt)
# STREAM loads rpy2 when it is imported and rpy2 has to know where R is installed
os.environ["R_HOME"] = subprocess.check_output(["R", "RHOME"], text=True).strip()
sys.path.insert(0, dataPath + "Software/STREAM-d20cc1faea58df10c53ee72447a9443f4b6c8e03/")
import stream


# the trajectories were computed with slingshot (01_Trajectories_Slingshot_Analysis.R)
# the plot function of STREAM is used to show the three lineages as a tree with the cells next to the branches
# STREAM expects its own tree object, so the slingshot result is put into this format first

clusters = ["OPC", "OL1", "OL2", "OL3", "OL4", "OL5", "OL6"]
cluster_color = {"OPC": "#F2CC5E", "OL1": "#6FA8DC", "OL2": "#E06666", "OL3": "#95D5B2", "OL4": "#B39BE8", "OL5": "#D9826B", "OL6": "#3B86C4"}

# read tables of the slingshot analysis
# cells: cell_id, sample_id, diagnosis, final_cluster_label, assigned_lineage, assigned_lineage_end_cluster, assigned_weight, assigned_pseudotime, PCA3_1, PCA3_2, PCA3_3
cells = pd.read_csv(tablePath + "Trajectories_Cells.txt.gz", sep="\t")
# curves: lineage, terminal_cluster, curve_order, curve_lambda, PCA3_1, PCA3_2, PCA3_3, ...
curves = pd.read_csv(tablePath + "Trajectories_Curves.txt", sep="\t")
# attached: position of each cell on its curve (attached_PCA3_1, attached_PCA3_2, attached_PCA3_3)
attached = pd.read_csv(tablePath + "Trajectories_CellsOnCurves.txt.gz", sep="\t")
# paths: lineage, path_order, cluster, root_cluster, terminal_cluster
paths = pd.read_csv(tablePath + "Trajectories_LineagePaths.txt", sep="\t")

lineages = {}
for name, frame in paths.groupby("lineage", sort=False):
    lineages[str(name)] = frame.sort_values("path_order")["cluster"].astype(str).tolist()


#
### Pseudotime of the nodes of the tree (one node per cluster)
#

# inner nodes: the curve point that is closest to the center of the cluster, mean over the lineages that pass the cluster
centers = cells.groupby("final_cluster_label", observed=True)[["PCA3_1", "PCA3_2", "PCA3_3"]].median()
inner_nodes = {"OL2": ["Lineage 1", "Lineage 2", "Lineage 3"], "OL4": ["Lineage 1", "Lineage 2"], "OL1": ["Lineage 1", "Lineage 2"]}
node_time = {"OPC": 0.0}
for node, node_lineages in inner_nodes.items():
    target = centers.loc[node].to_numpy(dtype=float)
    projected = []
    for lineage in node_lineages:
        curve = curves[curves["lineage"].astype(str) == lineage].sort_values("curve_order")
        xyz = curve[["PCA3_1", "PCA3_2", "PCA3_3"]].to_numpy(dtype=float)
        nearest = int(np.argmin(np.linalg.norm(xyz - target, axis=1)))
        projected.append(float(curve.iloc[nearest]["curve_lambda"]))
    node_time[node] = float(np.mean(projected))

# end nodes: end of the curve
for terminal, lineage in {"OL5": "Lineage 1", "OL6": "Lineage 2", "OL3": "Lineage 3"}.items():
    node_time[terminal] = float(curves.loc[curves["lineage"].astype(str) == lineage, "curve_lambda"].max())

# edges of the tree
edges = []
for path in lineages.values():
    for edge in zip(path[:-1], path[1:]):
        if edge not in edges:
            edges.append(edge)


#
### Assign each cell to an edge of the tree
#

adapted = cells.merge(attached[["cell_id", "attached_PCA3_1", "attached_PCA3_2", "attached_PCA3_3"]], on="cell_id", how="inner", validate="one_to_one")

# distance of the cell to its curve
original_xyz = adapted[["PCA3_1", "PCA3_2", "PCA3_3"]].to_numpy(dtype=float)
attached_xyz = adapted[["attached_PCA3_1", "attached_PCA3_2", "attached_PCA3_3"]].to_numpy(dtype=float)
adapted["branch_dist"] = np.linalg.norm(original_xyz - attached_xyz, axis=1)

# edge: by the pseudotime of the cell along the path of its lineage
adapted["graph_edge"] = ""
for lineage, path in lineages.items():
    indices = adapted.index[adapted["assigned_lineage"].astype(str) == lineage].to_numpy()
    end_time = np.asarray([node_time[node] for node in path[1:]], dtype=float)
    pseudotime = adapted.loc[indices, "assigned_pseudotime"].to_numpy(dtype=float)
    edge_index = np.searchsorted(end_time, pseudotime, side="left")
    edge_index = np.clip(edge_index, 0, len(path) - 2)
    adapted.loc[indices, "graph_edge"] = [f"{path[i]}->{path[i + 1]}" for i in edge_index]

# cell annotation in the format of STREAM
edge_lookup = {f"{start}->{end}": (start, end) for start, end in edges}
starts = adapted["graph_edge"].str.split("->").str[0]

obs = pd.DataFrame(index=adapted["cell_id"].astype(str))
obs.index.name = "cell_id"
obs["final_cluster_label"] = pd.Categorical(adapted["final_cluster_label"].astype(str).to_numpy(), categories=clusters, ordered=True)
obs["branch_id"] = pd.Series(adapted["graph_edge"].map(edge_lookup).to_list(), index=obs.index, dtype=object)
obs["branch_lam"] = adapted["assigned_pseudotime"].to_numpy(dtype=float) - starts.map(node_time).to_numpy(dtype=float)
obs["branch_dist"] = adapted["branch_dist"].to_numpy(dtype=float)
obs["OPC_pseudotime"] = adapted["assigned_pseudotime"].to_numpy(dtype=float)

# tree in the format of STREAM
flat_tree = nx.Graph()
for node in node_time:
    flat_tree.add_node(node, label=node)
for start, end in edges:
    flat_tree.add_edge(start, end, id=(start, end), len=float(node_time[end] - node_time[start]), color="#25211F")
# STREAM still uses the old name of the node dictionary of networkx
flat_tree.node = flat_tree._node

adata = ad.AnnData(X=np.empty((len(obs), 0), dtype=np.float32), obs=obs)
adata.uns["flat_tree"] = flat_tree
adata.uns["final_cluster_label_color"] = cluster_color


#
### 4500 cells are shown, sampled proportional to the number of cells of each edge
#

n_show = 4500
counts = adapted["graph_edge"].value_counts(sort=False).sort_index()
raw = counts / counts.sum() * n_show
n_per_edge = np.floor(raw).astype(int)
n_per_edge[n_per_edge < 1] = 1
# the cells that are missing to 4500 because of rounding go to the edges with the largest remainder
remainder = n_show - int(n_per_edge.sum())
for edge in (raw - np.floor(raw)).sort_values(ascending=False).index[:remainder]:
    n_per_edge.loc[edge] += 1

pieces = []
for index, (edge, n_take) in enumerate(n_per_edge.items()):
    subset = adapted[adapted["graph_edge"] == edge]
    pieces.append(subset.sample(n=int(n_take), replace=False, random_state=100 + index))
selected = pd.concat(pieces, ignore_index=True)
selected = selected.sample(frac=1.0, random_state=100).reset_index(drop=True)

display = adata[selected["cell_id"].astype(str).tolist()].copy()

# positions of the cells in the plot
stream.add_stream_sc_pos(display, root="OPC", dist_scale=0.3, dist_pctl=95, preference=["OL5", "OL6", "OL3"])

coordinates = pd.DataFrame(display.obsm["X_stream_OPC"], columns=["pseudotime", "vertical_position"])
coordinates.insert(0, "cell_id", display.obs_names.astype(str))
coordinates.to_csv(tablePath + "Trajectories_STREAM_CellPositions.txt.gz", sep="\t", index=False)
selected[["cell_id", "graph_edge", "final_cluster_label"]].to_csv(tablePath + "Trajectories_STREAM_ShownCells.txt", sep="\t", index=False)

# the cells that are shown, with the tree, for the plot (03_Trajectories_STREAM_Plots.py)
with open(objectPath + "Trajectories_STREAM_PlotData.pkl", "wb") as handle:
    pickle.dump(display, handle)
