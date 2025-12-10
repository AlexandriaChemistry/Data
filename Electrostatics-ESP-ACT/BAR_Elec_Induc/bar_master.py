#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

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
    width = 0.12

    plt.figure(figsize=(18, 8), dpi=300)

    plt.bar(x - 3*width, SAPT, width, label="SAPT2+(CCD)-δMP2", color="#5e3c99", edgecolor="black", linewidth=5)
    plt.bar(x - 2*width, TIP4P, width, label="TIP4P-Ew",       color="#c8f7c5", edgecolor="black", linewidth=0.7)
    plt.bar(x - width,   MBIS, width, label="MBIS-S",          color="#ffb000", edgecolor="black", linewidth=0.7)
    plt.bar(x,           CHARMM, width, label="CHARMM Drude",  color="#c9b7ff", edgecolor="black", linewidth=0.7)
    plt.bar(x + width,   GV4, width, label="PC+GV4",           color="#e41a1c", edgecolor="black", linewidth=5)
    plt.bar(x + 2*width, GS4, width, label="PC+GS4",           color="#00a9cf", edgecolor="black", linewidth=5)

    plt.xticks(x, ions_display, fontsize=30, rotation=20)
    plt.yticks(fontsize=30)
    plt.ylabel("Electrostatic Energy (kJ/mol)", fontsize=30)
    plt.axhline(0, color="black", linewidth=1.2)
    plt.grid(axis="y", linestyle="--", alpha=0.25)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.legend(fontsize=20, frameon=False, ncol=3, bbox_to_anchor=(0.5, -0.22), loc="upper center")
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.05, left=0.05)
    plt.savefig("table4.pdf", dpi=600, bbox_inches="tight")
    plt.close()


def plot_table6(texfile):
    ions, data = readLATEX(texfile)
    r, SAPT, PC, MBIS, Walz, GV4, GS4 = data.T
    ions_display = ["Li⁺–F⁻", "Li⁺–Cl⁻", "Li⁺–Br⁻", "Na⁺–F⁻", "Na⁺–Cl⁻", "Na⁺–Br⁻", "K⁺–F⁻", "K⁺–Cl⁻", "K⁺–Br⁻"]
    x = np.arange(len(ions_display))
    width = 0.13

    plt.figure(figsize=(18, 7), dpi=300)

    plt.bar(x - 2.5*width, SAPT, width, label="SAPT2+(CCD)-δMP2", color="#5e3c99", edgecolor="black", linewidth=5)
    plt.bar(x - 1.5*width, PC,   width, label="PC",               color="#c9b7ff", edgecolor="black", linewidth=1)
    plt.bar(x - 0.5*width, MBIS, width, label="MBIS-S",           color="#ffb000", edgecolor="black", linewidth=1)
    plt.bar(x + 0.5*width, Walz, width, label="Walz Model",       color="#c6f9e3", edgecolor="black", linewidth=1)
    plt.bar(x + 1.5*width, GV4,  width, label="PC+GV4",           color="#e31a1c", edgecolor="black", linewidth=5)
    plt.bar(x + 2.5*width, GS4,  width, label="PC+GS4",           color="#00a9cf", edgecolor="black", linewidth=5)

    plt.xticks(x, ions_display, fontsize=30, rotation=20)
    plt.yticks(fontsize=30)
    plt.ylabel("Electrostatic Energy (kJ/mol)", fontsize=30)
    plt.axhline(0, color="black", linewidth=1.2)
    plt.grid(axis="y", linestyle="--", alpha=0.25)

    ax = plt.gca()
    ax.spines["bottom"].set_linewidth(1.8)
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.legend(fontsize=20, frameon=False, ncol=3, bbox_to_anchor=(0.5, -0.25), loc="upper center")
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.05, left=0.05)
    plt.savefig("table6.pdf", dpi=1200, bbox_inches="tight")
    plt.close()


def plot_table5(texfile):
    ions, data = readLATEX(texfile)
    r, SAPT, CHARMM, GS4 = data.T

    ions_display = ["Li⁺–H₂O", "Na⁺–H₂O", "K⁺–H₂O", "F⁻–H₂O", "Cl⁻–H₂O", "Br⁻–H₂O", "H₂O–H₂O"]
    x = np.arange(len(ions_display))
    width = 0.24

    plt.figure(figsize=(14, 6), dpi=300)

    plt.bar(x - width, SAPT, width, label="SAPT2+(CCD)-δMP2", color="#5e3c99", edgecolor="black", linewidth=5)
    plt.bar(x,        CHARMM, width, label="CHARMM Drude",     color="#e31a1c", edgecolor="black", linewidth=1)
    plt.bar(x + width, GS4,   width, label="PC+GS4",            color="royalblue", edgecolor="black", linewidth=5)

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

    plt.legend(fontsize=20, frameon=False)
    plt.tight_layout()
    plt.savefig("table5.pdf", dpi=1200, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    plot_table4("ion-water-SAPT2-TIP4Pew-ACT4S.tex")
    plot_table6("Ions-sapt2-JC-Walz2018a-ACT.tex")
    plot_table5("ion-water-induction.tex")
    print("figs generated: table4.pdf, table6.pdf, table5.pdf")

