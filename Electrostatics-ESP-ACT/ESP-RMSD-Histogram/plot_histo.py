#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

x, y = [], []
filename = "histo-HF-6-311G**.xvg"
with open(filename) as f:
    for line in f:
        if line.startswith("@") or line.strip() == "":
            continue
        parts = line.split()
        x.append(float(parts[0]))
        y.append(float(parts[1]))
x = np.array(x)
y = np.array(y)

fig, ax = plt.subplots(figsize=(10,6))
bar_width = (x[1]-x[0])*0.8

for i in range(len(x)):
    ax.bar(x[i], y[i], width=bar_width, color="royalblue",
           edgecolor="cornflowerblue", linewidth=1.5, alpha=0.95, zorder=3)

avg_rmsd = 5.7
#ax.axvline(avg_rmsd, color="#DC143C", lw=4, linestyle='--', zorder=6)
#ax.axvspan(avg_rmsd-0.05, avg_rmsd+0.05, color="#DC143C", alpha=0.25, zorder=6)
#ax.text(avg_rmsd+0.05, max(y)*0.92, f"Average RMSD = {avg_rmsd:.1f} kJ/mol e",
#        color="#DC143C", fontsize=22, fontweight='bold', rotation=90, va='top', zorder=7)

ax.set_xlabel("RMSD (kJ/mol e)", fontsize=32)
ax.set_ylabel("Frequency (a.u.)", fontsize=32)

ax.tick_params(axis='both', which='major', labelsize=32, width=2, length=7)
ax.tick_params(axis='both', which='minor', labelsize=32, width=1.5, length=4)

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1.5)

ax.set_xlim(min(x)-0.05, max(x)+0.05)
ax.set_ylim(0, max(y)*1.35)

plt.tight_layout()
plt.savefig("rmsd_histo.pdf", dpi=300)
plt.show()

