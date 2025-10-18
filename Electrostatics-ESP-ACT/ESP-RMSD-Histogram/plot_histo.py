#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from matplotlib.colors import LinearSegmentedColormap

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
max_idx = np.argmax(y)

for i in range(len(x)):
    color = "#DC143C" if abs(i - max_idx) <= 1 else "#4169E1"
    ax.bar(x[i], y[i], width=bar_width, color=color, edgecolor="black", linewidth=1.5, alpha=0.95, zorder=3)

x_smooth = np.linspace(x.min(), x.max(), 2000)
spline = make_interp_spline(x, y, k=5)
y_smooth = spline(x_smooth)
ax.fill_between(x_smooth, 0, y_smooth, color='grey', alpha=0.1, zorder=1)
ax.plot(x_smooth, y_smooth, color='black', lw=2.5, alpha=0.9, zorder=4)

ax.scatter(x, y, color='black', s=60, zorder=5)

avg_rmsd = 5.7
ax.axvline(avg_rmsd, color="#DC143C", lw=4, linestyle='--', zorder=6, label=f"Average RMSD = {avg_rmsd:.1f} kJ/mol e")
ax.axvspan(avg_rmsd-0.05, avg_rmsd+0.05, color="#DC143C", alpha=0.25, zorder=6)
#ax.text(avg_rmsd+0.02, max(y)*0.9 color="#DC143C",
#        fontsize=22, fontweight='bold', rotation=90, va='top', zorder=7)

ax.set_xlabel("RMSD (kJ/mol e)", fontsize=32)
ax.set_ylabel("Frequency (a.u.)", fontsize=32)

ax.tick_params(axis='both', which='major', labelsize=32, width=2, length=7)
ax.tick_params(axis='both', which='minor', labelsize=32, width=1.5, length=4)

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1.5)

ax.set_xlim(min(x)-0.05, max(x)+0.05)
ax.set_ylim(0, max(y)*1.35)

plt.legend(fontsize=28)
plt.tight_layout()
plt.savefig("rmsd_histo.pdf", dpi=300)
plt.show()

