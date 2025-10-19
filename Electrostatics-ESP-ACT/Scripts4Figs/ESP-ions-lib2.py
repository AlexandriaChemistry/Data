#!/usr/bin/env python3
import os, sys, json
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../Scripts4Tabs')
from walz2018a_SAPT0_P1S2S import one_4pi_eps0

input_json = 'output_4_100.json'
with open(input_json, 'r') as f:
    output_data = json.load(f)

start = 80
delta = 5
xref = output_data["data"]["Br"][0]
indices = list(range(start, len(xref), delta))
distances = np.array([xref[i] for i in indices])

pc = [0]*start
for i in range(start, len(xref)):
    pc.append(one_4pi_eps0 / xref[i] if abs(xref[i]) > 1e-8 else 0.0)
pc = np.array(pc)

plt.rcParams.update({
    "axes.linewidth": 2,
    "axes.edgecolor": "black",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
})

ion_colors = {
    "F": "#0077be",
    "Cl": "#2ca02c",
    "Br": "blue",
    "Li": "#ff7f0e",
    "Na": "#f1c40f",
    "K": "#d62728",
    "PC+": "black",
    "PC-": "black"
}

ion_markers = {
    "F": "o",
    "Cl": "s",
    "Br": "D",
    "Li": "^",
    "Na": "v",
    "K": "*",
    "PC+": ".",
    "PC-": "1"
}

marker_size = 360
edge_width = 2.2

fig, ax = plt.subplots(figsize=(12, 8))

for ion in ["F", "Cl", "Br"]:
    yvals = np.array([output_data["data"][ion][1][i] for i in indices])
    ax.scatter(
        distances, yvals,
        color=ion_colors[ion], s=marker_size,
        marker=ion_markers[ion], linewidth=edge_width, 
        alpha=0.9, label=f"{ion}-"
    )
ax.scatter(distances, -pc[indices], color='black', s=marker_size*0.9,
           marker=ion_markers["PC-"], linewidth=edge_width+0.5, label=f"PC-")

for ion in ["Li", "Na", "K"]:
    yvals = np.array([output_data["data"][ion][1][i] for i in indices])
    ax.scatter(
        distances, yvals,
        color=ion_colors[ion], s=marker_size,
        marker=ion_markers[ion], linewidth=edge_width,
        alpha=0.9, label=f"{ion}+"
    )
ax.scatter(distances, pc[indices], color='black', s=marker_size*0.9,
           marker=ion_markers["PC+"], linewidth=edge_width+0.5, label="PC+")

ax.tick_params(axis='both', which='major', labelsize=26,
               width=2, length=8, direction='inout')
ax.set_facecolor((0.98, 0.98, 0.98))
for spine in ax.spines.values():
    spine.set_linewidth(2)

ax.set_xlabel("Distance ($\\mathrm{\\AA}$)", fontsize=32, labelpad=14)
ax.set_ylabel("ESP (kJ/mol·e)", fontsize=32, labelpad=14)
ax.legend(fontsize=20, frameon=True, loc='upper right', edgecolor='black', fancybox=True)

os.makedirs("../Figures", exist_ok=True)
plt.tight_layout()
plt.savefig("../Figures/ion-esp.pdf", dpi=600, bbox_inches='tight')
plt.show()

