Data shown in the figures of the article: one folder per figure, tab-separated files named after the panel (for example Figure_1C_Volcano.txt.gz).
Each file holds only what is drawn in the panel (positions, colours, sizes, written values and labels), taken from the results of the scripts of the figure folder.
The files *_Nuclei (Figure 5D, 5E, Supplementary Figure 1B to 1E) hold the nuclei shown as points, without cell or donor identifiers.

Column names
- z_01 to z_50 (Figure 4B to 4D): 50 groups of HC nuclei with the same number of nuclei, ordered by pseudotime (Figure_4*_Columns.txt: position on the axis, number of nuclei and pseudotime of each group)
- HC_01 to HC_50 and disease_01 to disease_50 (Figure 5): 50 groups with the same number of nuclei (HC and disease together), ordered along OL pseudotime, for HC and for the disease (Figure_5A_Modules_Points.txt: position on the axis and median pseudotime of each group)
- axis_position: position on the horizontal axis of the panel, from 0 to 1
- empty fields: nothing is drawn there (for example no dot for a transcription factor that is not significant)

Not included are the UMAPs of Figure 3A, 3B and 3E, the nuclei of Figure 4A and the donor proportions of Supplementary Figure 1A
(single nuclei and donors of the snRNA-seq data set of Rydbirk et al., 2026).
Figure 1A, 1B, 2D and 2E are images and measurements of the HOG cells, not analyses of this repository.
