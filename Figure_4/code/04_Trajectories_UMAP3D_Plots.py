import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
import pandas as pd

# change to your project directory
home_directory = "path/to/project_directory/"

# output folder
outputPath = home_directory + "Figure_4/"
figurePath = outputPath + "figures/"
tablePath = outputPath + "tables/Trajectories/"

mpl.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "axes.unicode_minus": False})


clusters = ["OPC", "OL1", "OL2", "OL3", "OL4", "OL5", "OL6"]
cluster_color = {"OPC": "#F2CC5E", "OL1": "#6FA8DC", "OL2": "#E06666", "OL3": "#95D5B2", "OL4": "#B39BE8", "OL5": "#D9826B", "OL6": "#3B86C4"}
lineage_color = {"Lineage 1": "#D9826B", "Lineage 2": "#3B86C4", "Lineage 3": "#2E8B57"}

# sizes in points: all text is 10 pt, line widths and point sizes follow the size of the figure
textSize = 10
figureScale = 0.78

# read the three-dimensional UMAP of the cells and the slingshot curves in this UMAP (made by 01_Trajectories_Slingshot_Analysis.R)
# cells: cell_id, sample_id, diagnosis, final_cluster_label, UMAP3D_1, UMAP3D_2, UMAP3D_3
cells = pd.read_csv(tablePath + "Trajectories_UMAP3D_Cells.txt.gz", sep="\t", low_memory=False)
# curves: UMAP3D_1, UMAP3D_2, UMAP3D_3, segment_length, lineage, curve_order
curves = pd.read_csv(tablePath + "Trajectories_UMAP3D_Curves.txt", sep="\t")

cells["final_cluster_label"] = pd.Categorical(cells["final_cluster_label"], clusters, ordered=True)
curves = curves.sort_values(["lineage", "curve_order"]).reset_index(drop=True)

# range of the axes: all cells and curves with a small margin
columns = ["UMAP3D_1", "UMAP3D_2", "UMAP3D_3"]
combined = pd.concat([cells[columns], curves[columns]], ignore_index=True)
low = combined.min().to_numpy(dtype=float)
high = combined.max().to_numpy(dtype=float)
span = np.maximum(high - low, 1e-8)
low, high = low - span * 0.035, high + span * 0.035
span = high - low

fig = plt.figure(figsize=(6.35 * figureScale, 4.62 * figureScale), dpi=150)
ax = fig.add_axes([0.025, 0.135, 0.95, 0.84], projection="3d")
# the plotting area has no white background
ax.set_facecolor((1, 1, 1, 0))

# cells
for cluster in clusters:
    subset = cells[cells["final_cluster_label"].astype(str) == cluster]
    ax.scatter(subset["UMAP3D_1"], subset["UMAP3D_2"], subset["UMAP3D_3"], c=cluster_color[cluster], s=0.72 * figureScale ** 2, alpha=0.25,
               depthshade=True, linewidths=0, rasterized=True, zorder=2)

# curves of the three lineages
curveCollections = []
for lineage, subset in curves.groupby("lineage", sort=False):
    subset = subset.sort_values("curve_order")
    xyz = subset[columns].to_numpy(dtype=float)
    segments = np.stack([xyz[:-1], xyz[1:]], axis=1)
    collection = Line3DCollection(segments, colors=lineage_color[str(lineage)], linewidths=2.20 * figureScale, alpha=0.98, linestyles="solid")
    collection.set_zorder(8)
    ax.add_collection3d(collection)
    curveCollections.append(collection)

# cluster names at the median position of the cells, shifted a bit so that the names are readable
label_shift = {"OPC": (0.02, 0.03, 0.08), "OL1": (-0.06, -0.01, 0.08), "OL2": (0.04, 0.06, 0.06), "OL3": (0.03, -0.06, 0.10),
               "OL4": (0.03, -0.04, -0.09), "OL5": (-0.06, -0.02, 0.06), "OL6": (-0.02, -0.05, -0.10)}
centers = cells.groupby("final_cluster_label", observed=True)[columns].median()
for cluster in clusters:
    xyz = centers.loc[cluster].to_numpy(dtype=float) + np.asarray(label_shift[cluster], dtype=float) * span
    ax.text(xyz[0], xyz[1], xyz[2], cluster, fontsize=textSize, fontweight="normal", color="black", zorder=12)

# axes and view
ax.set_xlim(low[0], high[0])
ax.set_ylim(low[1], high[1])
ax.set_zlim(low[2], high[2])
ax.set_box_aspect(tuple(span))
ax.view_init(elev=22.0, azim=-60.0, roll=0.0)
ax.set_proj_type("persp", focal_length=0.72)

ax.set_xlabel("UMAP1", labelpad=-12, fontsize=textSize, fontweight="normal", color="black")
ax.set_ylabel("UMAP2", labelpad=-12, fontsize=textSize, fontweight="normal", color="black")
ax.set_zlabel("UMAP3", labelpad=-12, fontsize=textSize, fontweight="normal", color="black")
ax.set_xticks(np.linspace(low[0], high[0], 4))
ax.set_yticks(np.linspace(low[1], high[1], 4))
ax.set_zticks(np.linspace(low[2], high[2], 4))
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.tick_params(axis="both", which="major", length=0, pad=0)
for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
    axis.pane.set_facecolor((1, 1, 1, 0))
    axis.pane.set_edgecolor("#D8D2CB")
    axis.line.set_color("#4A4642")
    axis.line.set_linewidth(1.0 * figureScale)
    axis._axinfo["grid"]["color"] = (0.84, 0.81, 0.77, 0.70)
    axis._axinfo["grid"]["linewidth"] = 0.55 * figureScale
    axis._axinfo["tick"]["inward_factor"] = 0.0
    axis._axinfo["tick"]["outward_factor"] = 0.0

# the cells keep the order of mplot3d (sorted by depth), the curves and the cluster names are drawn over all cells
fig.canvas.draw()
ax.computed_zorder = False
cellLayer = max(collection.get_zorder() for collection in ax.collections if collection not in curveCollections)
for k, collection in enumerate(curveCollections):
    collection.set_zorder(cellLayer + 1 + k)
for artist in ax.texts:
    artist.set_zorder(cellLayer + 10)
# legend lines 0.5 pt and 18 pt long, as in the legends of Figure 5
handles = [Line2D([0], [0], color=color, lw=0.5, solid_capstyle="round", label=lineage) for lineage, color in lineage_color.items()]
legend = fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.50, 0.075), ncol=3, frameon=False, handlelength=1.8, handletextpad=0.55,
           columnspacing=1.2, prop={"family": "Arial", "size": textSize, "weight": "normal"})

# the files are cropped to what is drawn: panes, axis lines, axis titles and the legend (mplot3d would keep the empty margin of its plotting area);
# cells, curves and cluster names lie inside the panes
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
drawn = [axis.pane for axis in (ax.xaxis, ax.yaxis, ax.zaxis)] + [axis.line for axis in (ax.xaxis, ax.yaxis, ax.zaxis)] + [axis.label for axis in (ax.xaxis, ax.yaxis, ax.zaxis)]
drawn += [legend]
crop = Bbox.union([artist.get_window_extent(renderer) for artist in drawn]).transformed(fig.dpi_scale_trans.inverted()).padded(0.04)
fig.savefig(figurePath + "Figure_4A_Trajectories_UMAP3D.pdf", facecolor="white", transparent=False, bbox_inches=crop)

# the same plot with the cells as vector graphic
for collection in ax.collections:
    collection.set_rasterized(False)
fig.savefig(figurePath + "Figure_4A_Trajectories_UMAP3D_FullVector.pdf", facecolor="white", transparent=False, bbox_inches=crop)
plt.close(fig)
