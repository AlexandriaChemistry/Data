#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

lw    = 0.3
legfs = 28
#colors = [ "#5e3c99", "#c8f7c5", "#ffb000", "#0000cc", "#e41a1c", "#e41a1c" ]
colors  = [ "#aaaaaa", "#00cc00",  "#00eeee", "#0000cc", "#ffff00", "#ff0000" ]
hatches = [ '+', '', '-', '', 'x', '', 'O', '.', '*', '\\' ]

def readLATEX(filepath):
    labels, rows = [], []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if "&" not in line or line.startswith("\\hline") or "RMSD" in line or "MSE" in line:
                continue
            parts = [p.strip() for p in line.replace("\\\\", "").split("&")]
            nums = []
            for entry in parts[1:]:
                clean = entry.replace("\\", "").replace(",", "")
                try:
                    nums.append(float(clean))
                except:
                    continue
            if nums:
                labels.append(parts[0])
                rows.append(nums)
    return labels, np.array(rows)


def plot_table4(texfile):
    ions, data = readLATEX(texfile)
    r, SAPT, TIP4P, MBIS, CHARMM, GV4, GS4 = data.T
    ions_display = ["Li⁺–H₂O", "Na⁺–H₂O", "K⁺–H₂O", "F⁻–H₂O", "Cl⁻–H₂O", "Br⁻–H₂O", "H₂O–H₂O"]
    x = np.arange(len(ions_display))
    mydata = [ { "data": SAPT,  "label": "SAPT2+(CCD)-δMP2" },
               { "data": TIP4P,  "label": "PC/TIP4P-Ew" },
               { "data": MBIS,   "label": "MBIS-S" },
               { "data": CHARMM, "label": "CHARMM Drude" },
               { "data": GV4,    "label": "PC+GV4" },
               { "data": GS4,    "label": "PC+GS4" } ]
    
    width = 0.12

    plt.figure(figsize=(12, 15), dpi=300)

    x0 = x - 3*width
    i  = 0
    for md in mydata:
        plt.barh(x0, md["data"],  width, label=md["label"], color=colors[i], edgecolor="black", linewidth=lw, hatch=hatches[i])
        x0 += width
        i += 1

    plt.yticks(x, ions_display, fontsize=30, rotation=20)
    plt.xticks(fontsize=20)
    plt.xlabel("Electrostatic Energy (kJ/mol)", fontsize=30)
#    plt.axhline(0, color="black", linewidth=1.2)
    plt.grid(axis="x", linestyle="--", alpha=0.25)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.legend(fontsize=legfs, frameon=False, ncol=2, bbox_to_anchor=(0.5, -0.07), loc="upper center")
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.05, left=0.05)
    plt.savefig("elec-water-ions.pdf", dpi=600, bbox_inches="tight")
    plt.close()


def plot_table6(texfile):
    ions, data = readLATEX(texfile)
    r, SAPT, PC, MBIS, Walz, GV4, GS4 = data.T
    ions_display = ["Li⁺–F⁻", "Li⁺–Cl⁻", "Li⁺–Br⁻", "Na⁺–F⁻", "Na⁺–Cl⁻", "Na⁺–Br⁻", "K⁺–F⁻", "K⁺–Cl⁻", "K⁺–Br⁻"]
    x = np.arange(len(ions_display))
    width = 0.13
    mydata = [ { "data": SAPT,   "label": "SAPT2+(CCD)-δMP2" },
               { "data": PC,     "label": "PC/TIP4P-Ew" },
               { "data": MBIS,   "label": "MBIS-S" },
               { "data": Walz,   "label": "Walz2018" },
               { "data": GV4,    "label": "PC+GV4" },
               { "data": GS4,    "label": "PC+GS4" } ]
    
    plt.figure(figsize=(12, 15), dpi=300)

    x0 = x - 2.5*width
    i  = 0
    for md in mydata:
        plt.barh(x0, md["data"],  width, label=md["label"], color=colors[i], edgecolor="black", linewidth=lw, hatch=hatches[i])
        x0 += width
        i += 1

    plt.yticks(x, ions_display, fontsize=30, rotation=20)
    plt.xticks(fontsize=20)
    plt.xlabel("Electrostatic Energy (kJ/mol)", fontsize=30)
#    plt.axhline(0, color="black", linewidth=1.2)
    plt.grid(axis="x", linestyle="--", alpha=0.25)

    ax = plt.gca()
    ax.spines["bottom"].set_linewidth(1.8)
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.legend(fontsize=legfs, frameon=False, ncol=2, bbox_to_anchor=(0.5, -0.07), loc="upper center")
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.05, left=0.05)
    plt.savefig("elec-ions.pdf", dpi=1200, bbox_inches="tight")
    plt.close()


def plot_table5(texfile):
    ions, data = readLATEX(texfile)
    r, SAPT, CHARMM, GS4 = data.T

    ions_display = ["Li⁺–H₂O", "Na⁺–H₂O", "K⁺–H₂O", "F⁻–H₂O", "Cl⁻–H₂O", "Br⁻–H₂O", "H₂O–H₂O"]
    x = np.arange(len(ions_display))
    width = 0.24

    plt.figure(figsize=(14, 6), dpi=300)

    plt.bar(x - width, SAPT, width, label="SAPT2+(CCD)-δMP2", color=colors[0], edgecolor="black", linewidth=lw, hatch=hatches[0])
    plt.bar(x,        CHARMM, width, label="CHARMM Drude",    color=colors[3], edgecolor="black", linewidth=lw, hatch=hatches[3])
    plt.bar(x + width, GS4,   width, label="PC+GS4",          color=colors[5], edgecolor="black", linewidth=lw, hatch=hatches[5])

    plt.xticks(x, ions_display, fontsize=30, rotation=20)
    plt.yticks(fontsize=30)
    plt.ylabel("Induction Energy (kJ/mol)", fontsize=28)
    plt.axhline(0, color="black", linewidth=1.2)
    plt.grid(axis="y", linestyle="--", alpha=0.25)

    ax = plt.gca()
    ax.spines["bottom"].set_linewidth(1.8)
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.legend(fontsize=legfs, frameon=False)
    plt.tight_layout()
    plt.savefig("induc-water-ions.pdf", dpi=1200, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    plot_table4("ion-water-SAPT2-TIP4Pew-ACT4S.tex")
    plot_table6("Ions-sapt2-JC-Walz2018a-ACT.tex")
    plot_table5("ion-water-induction.tex")

