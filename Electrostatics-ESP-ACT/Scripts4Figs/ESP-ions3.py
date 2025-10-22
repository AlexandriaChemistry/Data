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

# Larger text everywhere
plt.rcParams.update({
    "axes.linewidth": 2,
    "axes.edgecolor": "black",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "font.size": 28,              # Default text size
    "axes.titlesize": 34,         # Axis title size
    "axes.labelsize": 34,         # Label size
    "xtick.labelsize": 28,        # Tick label size
    "ytick.labelsize": 28,        # Tick label size
    "legend.fontsize": 28,        # Legend font size
})

ion_colors = {
    "F": "crimson",
    "Cl": "mediumseagreen",
    "Br": "cornflowerblue",
    "Li": "crimson",
    "Na": "mediumseagreen",
    "K": "cornflowerblue",
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
    "PC-": "."
}

marker_size = 360
edge_width = 2.2

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 14), sharex=False)

for ion in ["F", "Cl", "Br"]:
    yvals = np.array([output_data["data"][ion][1][i] for i in indices])
    ax1.scatter(
        distances, yvals,
        color=ion_colors[ion], s=marker_size,
        marker=ion_markers[ion],
        linewidth=edge_width, alpha=0.9, label=f"{ion}-"
    )
ax1.scatter(distances, -pc[indices], color='black', s=marker_size*1.2,
            marker=ion_markers["PC-"], linewidth=edge_width+0.5, label=f"PC-")

for ion in ["Li", "Na", "K"]:
    yvals = np.array([output_data["data"][ion][1][i] for i in indices])
    ax2.scatter(
        distances, yvals,
        color=ion_colors[ion], s=marker_size,
        marker=ion_markers[ion],
        linewidth=edge_width, alpha=0.9, label=f"{ion}+"
    )
ax2.scatter(distances, pc[indices], color='black', s=marker_size*1.2,
            marker=ion_markers["PC+"], linewidth=edge_width+0.5, label="PC+")

for ax, label in zip([ax1, ax2], ["A", "B"]):
    ax.text(-0.08, 1.1, label, transform=ax.transAxes,
            fontsize=38, va='top', ha='right', fontweight='bold')
    ax.tick_params(axis='both', which='major', width=2, length=8, direction='inout')
    ax.set_facecolor((0.98, 0.98, 0.98))
    ax.set_xlim(0.7, 3.)  # Limit x-axis to 4 Å
    for spine in ax.spines.values():
        spine.set_linewidth(2)

ax1.set_ylabel("ESP (kJ/mol·e)", fontsize=42, labelpad=18)
ax2.set_ylabel("ESP (kJ/mol·e)", fontsize=42, labelpad=18)
ax2.set_xlabel("Distance ($\\mathrm{\\AA}$)", fontsize=42, labelpad=18)

ax1.legend(frameon=True, loc='upper right', edgecolor='black', fontsize=28)
ax2.legend(frameon=True, loc='upper right', edgecolor='black', fontsize=28)
ax1.tick_params(axis='both', which='major', labelsize=42)
ax2.tick_params(axis='both', which='major', labelsize=42)

os.makedirs("../Figures", exist_ok=True)
plt.subplots_adjust(hspace=0., wspace=0.)
plt.tight_layout(h_pad=2.5)
plt.savefig("../Figures/ion-esp.pdf", dpi=600, bbox_inches='tight')
plt.show()

