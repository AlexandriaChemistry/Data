#!/usr/bin/env python
import numpy as np
import json

with open('../Charge_Models/sapt2_data.json', 'r') as file:
    data = json.load(file)

allowed_pairs = {
    "water#lithium",  
    "water#sodium",   
    "water#potassium",
    "water#fluoride", 
    "water#chloride", 
    "water#bromide",  
    "water#water"
}

data = {k: v for k, v in data.items() if k in allowed_pairs}

ions = [k.split('#')[1] if k != "water#water" else "water" for k in data]
rmins = [data[k]["rmin"] for k in data]
sapt_vals = [data[k]["eelec-SAPT"] for k in data]
jc_vals = [data[k]["JC"] for k in data]
swm4_vals = [data[k]["SWM4"] for k in data]
gcpgv_vals = [data[k]["GC+PGV"] for k in data]
pcgvs_vals = [data[k]["PC+GVS"] for k in data]

sapt_array = np.array(sapt_vals)
jc_array = np.array(jc_vals)
swm4_array = np.array(swm4_vals)
gcpgv_array = np.array(gcpgv_vals)
pcgvs_array = np.array(pcgvs_vals)

def compute_rmsd(pred, ref):
    return np.sqrt(np.mean((pred - ref) ** 2))

def compute_mse(pred, ref):
    return np.mean(pred - ref)

rmsd_values = [
    compute_rmsd(jc_array, sapt_array),
    compute_rmsd(swm4_array, sapt_array),
    compute_rmsd(gcpgv_array, sapt_array),
    compute_rmsd(pcgvs_array, sapt_array)
]

mse_values = [
    compute_mse(jc_array, sapt_array),
    compute_mse(swm4_array, sapt_array),
    compute_mse(gcpgv_array, sapt_array),
    compute_mse(pcgvs_array, sapt_array)
]

file_path = "ion-water-SAPT2-TIP4Pew-ACT4S.tex"
with open(file_path, "w") as file:
    file.write("\\begin{table}[ht]\n")
    file.write("\\centering\n")
    file.write("\\caption{\\textbf{Water-ion energies at their energy minimum.} Minimum energy distance (\\AA) between ions and water oxygen/hydrogen from Experiment (ref.~\\citenum{Heyrovska2006a}), and minimized water dimer (ref.~\\citenum{temelso2011benchmark}). Electrostatic energies are reported in kJ/mol from the SAPT2+(CCD)-$\\delta$MP2 method with an aug-cc-pVTZ basis set, TIP4P-Ew~\\cite{Horn2004a} with point charges representing ions, and SWM4-NDP~\\cite{Lamoureux2006a} with ions due to Yu {\\em et al.}~\\cite{Yu2010a}, point core+Gaussian vsite (GC+PGV), and point charge + Gaussian vsite and shell (PC+GVS) using ACT.}\n")
    file.write("\\label{tab:ion_water2}\n")
    file.write("\\begin{tabular}{lcccccc} \n")
    file.write("\\hline \n")
    file.write("Ion & r$_{min}$ & SAPT & TIP4P-Ew & SWM4-NDP & GC+PGV & PC+GVS\\\\\n")
    file.write("\\hline \n")
    for ion, rmin, sapt, jc, swm4, gcpgv, pcgvs in zip(ions, rmins, sapt_vals, jc_vals, swm4_vals, gcpgv_vals, pcgvs_vals):
        file.write(f"{ion} & {rmin:.2f} & {sapt:.1f} & {jc:.1f} & {swm4:.1f} & {gcpgv:.1f} & {pcgvs:.1f} \\\\\n")
    file.write("\\hline\n")
    file.write(f"RMSD & & & {rmsd_values[0]:.1f} & {rmsd_values[1]:.1f} & {rmsd_values[2]:.1f} & {rmsd_values[3]:.1f}\\\\\n")
    file.write(f"MSE & & & {mse_values[0]:.1f} & {mse_values[1]:.1f} & {mse_values[2]:.1f} & {mse_values[3]:.1f} \\\\\n")
    file.write("\\hline\n")
    file.write("\\end{tabular} \n")
    file.write("\\end{table}\n")

print(f"Please check file {file_path}")
