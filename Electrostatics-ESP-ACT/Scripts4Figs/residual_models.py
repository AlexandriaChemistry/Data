#!/usr/bin/env python3

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from potential_elec_functions import one_4pi_eps0, Point_core_gaussian_shell

FIG_DIR = "Figures"


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

    func_index_to_name = {0: "P+G"}
    func_index_to_function = {0: Point_core_gaussian_shell}

    fig, all_axes = plt.subplots(3, 3, figsize=(12, 10), sharex=False, sharey=False)
    axes = all_axes.flatten()
    idx = 0

    for row, cation in enumerate(cations):
        for col, anion in enumerate(anions):
            compound = cation + anion
            data_file = f'{data_dir}/distances_Electrostatics-{cation.lower()}{anion.lower()}.txt'

            dataset = np.loadtxt(data_file)
            x, y = dataset[:, 0], dataset[:, 1]

            sorted_idx = np.argsort(x)
            x_sorted = x[sorted_idx]
            y_sorted = y[sorted_idx]

            ener_pc = -one_4pi_eps0 / x_sorted
            ax = axes[idx]

            ax.text(
                0.5, 0.5, compound,
                transform=ax.transAxes,
                fontsize=20, ha='center'
            )

            ax.scatter(
                x_sorted, y_sorted - ener_pc,
                color="cornflowerblue", #colors[compound],
                label="PC",
                alpha=0.8,
                #edgecolor='k',
                s=120
            )

            for func_index in func_index_to_function:
                model_func = func_index_to_function[func_index]

                model_energies = model_func(
                    x_sorted,
                    params[cation][f"q_c_{func_index}"],
                    params[cation][f"q_s_{func_index}"],
                    params[anion][f"q_c_{func_index}"],
                    params[anion][f"q_s_{func_index}"],
                    params[cation][f"z2_{func_index}"],
                    params[anion][f"z2_{func_index}"]
                )

                diff = y_sorted - model_energies

                ax.plot(
                    x_sorted, diff,
                    label=func_index_to_name[func_index],
                    color="crimson", #"#7B3F61",
                    marker='^',
                    linestyle='None',
                    alpha=0.8,
                    markersize=12
                )

            if row < 2:
                ax.set_xticklabels([])
            if col != 0:
                ax.set_yticklabels([])

            ax.tick_params(axis='x', labelsize=24)
            ax.tick_params(axis='y', labelsize=24, pad=8)
            ax.legend(fontsize=14, loc='best', frameon=False)

            idx += 1

    first_row_axes = axes[0:3]
    y_min, y_max = np.inf, -np.inf
    for ax in first_row_axes:
        ylim = ax.get_ylim()
        y_min = min(y_min, ylim[0])
        y_max = max(y_max, ylim[1])
    for ax in first_row_axes:
        ax.set_ylim(y_min, y_max)

    fig.text(0.5, 0.04, 'Distance ($\\mathrm{\\AA}$)', ha='center', fontsize=24)
    fig.text(0.04, 0.5, 'Residual (kJ/mol)', va='center', rotation='vertical', fontsize=24)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0., wspace=0., left=0.2, bottom=0.12, right=0.95, top=0.95)

    output_file = "Residual-ESP.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    os.makedirs(FIG_DIR, exist_ok=True)
    for T in [10]:
        main(T)

