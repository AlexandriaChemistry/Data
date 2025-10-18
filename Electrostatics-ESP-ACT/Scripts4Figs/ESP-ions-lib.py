#!/usr/bin/env python3
import os, sys
import json
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
    if abs(xref[i]) > 1e-8:
        pc.append(one_4pi_eps0 / xref[i])
    else:
        pc.append(0.0)
pc = np.array(pc)

ion_colors = {
    "F": "dodgerblue",
    "Cl": "green",
    "Br": "royalblue",
    "Li": "orange",
    "Na": "gold",
    "K": "crimson",
    "PC+": "black",
    "PC-": "black"
}

marker_size = 250
edge_width = 1.5

plt.figure(figsize=(12, 7))

for ion in ["F", "Cl", "Br"]:
    yvals = np.array([output_data["data"][ion][1][i] for i in indices])
    plt.scatter(
        distances, yvals,
        color=ion_colors[ion], s=marker_size, edgecolor='black',
        linewidth=edge_width, label=f"{ion}-"
    )

plt.scatter(distances, -pc[indices], color='black', s=marker_size, marker='x', linewidth=edge_width, label="-PC")

for ion in ["Li", "Na", "K"]:
    yvals = np.array([output_data["data"][ion][1][i] for i in indices])
    plt.scatter(
        distances, yvals,
        color=ion_colors[ion], s=marker_size, edgecolor='black',
        linewidth=edge_width, label=f"{ion}+"
    )

plt.scatter(distances, pc[indices], color='black', s=marker_size, marker='x', linewidth=edge_width, label="PC+")

plt.xlabel("Distance ($\\mathrm{\\AA}$)", fontsize=32)
plt.ylabel("ESP (kJ/mol e)", fontsize=32)

plt.xticks(fontsize=32)
plt.yticks(fontsize=32)

plt.legend(fontsize=18, frameon=True, shadow=True, loc='upper right')

os.makedirs("../Figures", exist_ok=True)
plt.tight_layout()
plt.savefig("../Figures/ion-esp-absolute.pdf", dpi=300)
plt.show()

