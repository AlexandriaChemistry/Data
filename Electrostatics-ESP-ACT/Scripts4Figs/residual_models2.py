#!/usr/bin/env python3
import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.ndimage import gaussian_filter1d
from potential_elec_functions import one_4pi_eps0, Point_core_gaussian_shell

FIG_DIR = "Figures"

def glow_line(ax, x, y, color="royalblue", lw=3.5, glow_widths=[6, 10, 20], alpha_step=0.08):
    for i, w in enumerate(glow_widths[::-1]):
        alpha = alpha_step * (i + 1)
        ax.plot(x, y, color=color, linewidth=w, alpha=alpha, solid_capstyle='round', zorder=3)
    ax.plot(x, y, color=color, linewidth=lw, alpha=1.0, solid_capstyle='round', zorder=4)

def main(T: int):
    cations = ["Li", "Na", "K"]
    anions = ["F", "Cl", "Br"]
    data_dir = "SAPT_Alkali_Halides"

    with open(f'params_4_{T}.json', 'r') as json_f:
        params = json.load(json_f)

    colors = {
        "LiF": "crimson", "LiCl": "tomato", "LiBr": "lightcoral",
        "NaF": "forestgreen", "NaCl": "limegreen", "NaBr": "mediumseagreen",
        "KF": "royalblue", "KCl": "dodgerblue", "KBr": "cornflowerblue"
    }

    func_index_to_function = {0: Point_core_gaussian_shell}

    sns.set_theme(style="white")
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica"],
        "axes.linewidth": 1.5,
        "axes.labelsize": 28,
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "legend.fontsize": 14,
    })

    fig, all_axes = plt.subplots(3, 3, figsize=(12, 11))
    axes = all_axes.flatten()
    idx = 0

    first_row_ylims = []

    for row, cation in enumerate(cations):
        for col, anion in enumerate(anions):
            compound = cation + anion
            data_file = f'{data_dir}/distances_Electrostatics-{cation.lower()}{anion.lower()}.txt'
            dataset = np.loadtxt(data_file)
            x, y = dataset[:, 0], dataset[:, 1]
            sorted_idx = np.argsort(x)
            x_sorted, y_sorted = x[sorted_idx], y[sorted_idx]
            ener_pc = -one_4pi_eps0 / x_sorted
            ax = axes[idx]

            ax.text(0.5, 0.5, compound, transform=ax.transAxes,
                    fontsize=20, fontweight='bold', ha='center', va='center', color='black')

            model_func = func_index_to_function[0]
            model_energies = model_func(
                x_sorted,
                params[cation][f"q_c_0"],
                params[cation][f"q_s_0"],
                params[anion][f"q_c_0"],
                params[anion][f"q_s_0"],
                params[cation][f"z2_0"],
                params[anion][f"z2_0"]
            )
            diff = y_sorted - model_energies
            smooth = gaussian_filter1d(diff, sigma=1.2)
            glow_line(ax, x_sorted, smooth, color="black", lw=3, glow_widths=[6, 12, 20], alpha_step=0.05)
            glow_dummy = ax.plot([], [], color="black", lw=3, label="PG")[0]

            pg_scatter = ax.scatter(
                x_sorted, y_sorted - ener_pc,
                color=colors[compound],
                edgecolors='black',
                linewidths=0.6,
                s=120,
                alpha=0.85,
                label="PC"
            )

            ax.legend(handles=[glow_dummy, pg_scatter], loc='upper right', bbox_to_anchor=(1, 0.85), fontsize=16, frameon=False)

            if row == 0:
                first_row_ylims.append(ax.get_ylim())

            if row < 2:
                ax.set_xticklabels([])
            if col != 0:
                ax.set_yticklabels([])

            ax.tick_params(axis='both', which='major', labelsize=24)
            sns.despine(ax=ax, top=False, right=False)
            idx += 1

    if first_row_ylims:
        y_min = min([ylim[0] for ylim in first_row_ylims])
        y_max = max([ylim[1] for ylim in first_row_ylims])
        for ax in axes[0:3]:
            ax.set_ylim(y_min, y_max)

    fig.text(0.5, 0.04, 'Distance ($\\mathrm{\\AA}$)', ha='center', fontsize=28)
    fig.text(0.04, 0.5, 'Residual (kJ/mol)', va='center', rotation='vertical', fontsize=28)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0., wspace=0., left=0.15, bottom=0.12, right=0.98, top=0.95)
    output_file = f"{FIG_DIR}/Residual-ESP_AlignedLegend.pdf"
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    os.makedirs(FIG_DIR, exist_ok=True)
    for T in [10]:
        main(T)

