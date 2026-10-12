import os
import pickle
import subprocess
import sys
import xml.etree.ElementTree as ET

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrow, FancyArrowPatch
from matplotlib.text import Text

# change to your project directory
home_directory = "path/to/project_directory/"

# data folder
dataPath = home_directory + "Data/"

# output folder
outputPath = home_directory + "Figure_4/"
figurePath = outputPath + "figures/"
objectPath = outputPath + "objects/"

# STREAM 1.1 is not installed as a package, its source code is loaded from the data folder (download see Data/Software/README.txt)
# STREAM loads rpy2 when it is imported and rpy2 has to know where R is installed
os.environ["R_HOME"] = subprocess.check_output(["R", "RHOME"], text=True).strip()
sys.path.insert(0, dataPath + "Software/STREAM-d20cc1faea58df10c53ee72447a9443f4b6c8e03/")
import stream

clusters = ["OPC", "OL1", "OL2", "OL3", "OL4", "OL5", "OL6"]

# the cells that are shown, with the tree in the format of STREAM (made by 02_Trajectories_STREAM_Analysis.py)
with open(objectPath + "Trajectories_STREAM_PlotData.pkl", "rb") as handle:
    display = pickle.load(handle)


#
### Plot (Figure 4A)
#

# sizes in points: text is 10 pt (numbers at the axis 9 pt), line widths and point sizes follow the size of the figure
textSize = 10
tickTextSize = 9
lineWidth = 0.5
figureScale = 0.5

stream.set_figure_params(context="notebook", style="white", font="Arial", font_scale=1.25, dpi=100, dpi_save=300, figsize=[7 * figureScale, 4.5 * figureScale],
                         rc={"font.weight": "normal", "axes.labelweight": "normal", "axes.titleweight": "normal", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})

stream.plot_stream_sc(display, root="OPC", color=["final_cluster_label"], dist_scale=0.3, dist_pctl=95, preference=["OL5", "OL6", "OL3"],
                      fig_size=(7 * figureScale, 4.5 * figureScale), fig_legend_ncol=1, fig_legend_order={"final_cluster_label": clusters}, alpha=0.55, pad=1.08,
                      show_text=True, show_graph=True, save_fig=False, fig_path=figurePath, plotly=False)

fig = plt.gcf()
ax = fig.axes[0]

# the plot of STREAM is changed a bit:
# - cluster names of the nodes are placed by hand so that they do not cover the branches
# - no legend (the names are in the plot)
# - x-axis with ticks and an arrow
label_offset = {"OPC": {"offset": (5, 5), "ha": "left", "va": "bottom"}, "OL1": {"offset": (7, -7), "ha": "left", "va": "top"},
                "OL2": {"offset": (-7, -8), "ha": "right", "va": "top"}, "OL3": {"offset": (7, 4), "ha": "left", "va": "bottom"},
                "OL4": {"offset": (-6, 7), "ha": "right", "va": "bottom"}, "OL5": {"offset": (7, 4), "ha": "left", "va": "bottom"},
                "OL6": {"offset": (7, 5), "ha": "left", "va": "bottom"}}

node_position = {}
for artist in list(ax.texts):
    if artist.get_text() in clusters:
        node_position[artist.get_text()] = tuple(float(x) for x in artist.get_position())
        artist.remove()

ax.get_legend().remove()
ax.set_title("")
for collection in ax.collections:
    collection.set_sizes(np.asarray(collection.get_sizes()) * figureScale ** 2)
    collection.set_linewidth(np.asarray(collection.get_linewidth()) * figureScale)
for line in ax.lines:
    line.set_linewidth(line.get_linewidth() * figureScale)
    line.set_markersize(line.get_markersize() * figureScale)

for artist in list(ax.patches):
    if isinstance(artist, FancyArrow):
        artist.remove()
ax.spines["bottom"].set_visible(False)
ax.xaxis.set_ticks_position("bottom")
ax.tick_params(axis="x", which="major", bottom=True, top=False, labelbottom=True, direction="out", length=2, width=lineWidth, pad=2, color="black", labelcolor="black", labelsize=tickTextSize)
ax.add_patch(FancyArrowPatch((0.0, 0.0), (1.018, 0.0), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=11.5 * figureScale, linewidth=lineWidth, color="black", clip_on=False, zorder=20))

for label in clusters:
    ax.annotate(label, xy=node_position[label], xytext=label_offset[label]["offset"], textcoords="offset points", color="black",
                fontsize=textSize, fontfamily="Arial", fontweight="normal",
                ha=label_offset[label]["ha"], va=label_offset[label]["va"], annotation_clip=False)

# all text in Arial regular, black
for artist in fig.findobj(match=Text):
    if artist.get_text().strip():
        artist.set_fontfamily("Arial")
        artist.set_fontweight("normal")
        artist.set_color("black")
        if artist.get_text() not in clusters and not artist.get_text().lstrip("-").replace(".", "").isdigit():
            artist.set_fontsize(textSize)
ax.xaxis.label.set_fontsize(textSize)
for tick in ax.get_xticklabels():
    tick.set_fontsize(tickTextSize)
fig.canvas.draw()

fig.savefig(figurePath + "Figure_4A_Trajectories_STREAM.pdf", format="pdf", dpi=300, pad_inches=0.1, bbox_inches="tight")
plt.close(fig)
