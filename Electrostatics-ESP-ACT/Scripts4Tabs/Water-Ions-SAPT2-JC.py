#!/usr/bin/env python
import numpy as np
import json


with open('data.json', 'r') as file:
    data = json.load(file)


def rmsd(data, ref_key1, ref_key2):
    ref1_values = [values for values in data[ref_key1]]
    ref2_values = [values for values in data[ref_key2]]

    if len(ref1_values) != len(ref2_values):
        raise ValueError("Mismatched number of values between ref_key1 and ref_key2")

    ref1_values = np.array(ref1_values)
    ref2_values = np.array(ref2_values)

    return np.sqrt(np.mean((ref1_values - ref2_values)**2))



def mse(data, ref_key1, ref_key2):
    ref1_values = [values for values in data[ref_key1]]
    ref2_values = [values for values in data[ref_key2]]
    ref1_values = np.array(ref1_values)
    ref2_values = np.array(ref2_values)

    MSE = np.mean(np.subtract(ref1_values, ref2_values))
    return MSE




rmsd_value1 = rmsd(data, 'E$_{elec}$ (kJ/mol) 4W', 'E$_{elec}$ (kJ/mol) SAPT')
rmsd_value11 = rmsd(data, 'E$_{elec}$ (kJ/mol) SW', 'E$_{elec}$ (kJ/mol) SAPT')
rmsd_value2 = rmsd(data, 'E$_{elec}$ (kJ/mol) ACT$_{SS}$', 'E$_{elec}$ (kJ/mol) SAPT')
rmsd_value22 = rmsd(data, 'E$_{elec}$ (kJ/mol) ACT$_{GC}$', 'E$_{elec}$ (kJ/mol) SAPT')


mse_value1 = mse(data, 'E$_{elec}$ (kJ/mol) 4W', 'E$_{elec}$ (kJ/mol) SAPT')
mse_value11 = mse(data, 'E$_{elec}$ (kJ/mol) SW', 'E$_{elec}$ (kJ/mol) SAPT')
mse_value2 = mse(data, 'E$_{elec}$ (kJ/mol) ACT$_{SS}$', 'E$_{elec}$ (kJ/mol) SAPT')
mse_value22 = mse(data, 'E$_{elec}$ (kJ/mol) ACT$_{GC}$', 'E$_{elec}$ (kJ/mol) SAPT')

file_path = "ion-water-SAPT2-TIP4Pew-ACT4S.tex"

with open(file_path, "w") as file:
    file.write("\\begin{table}[ht]\n")
    file.write("\\centering\n")
    file.write("\\caption{\\textbf{Water-ion energies at their energy minimum.} Minimum energy distance (\\AA) between ions and water oxygen/hydrogen from Experiment (ref.~\\citenum{Heyrovska2006a}), and minimized water dimer (ref.~\\citenum{temelso2011benchmark}). Electrostatic energies are reported in kJ/mol from the SAPT2+(CCD)-$\\delta$MP2 method with an aug-cc-pVTZ basis set, TIP4P-Ew~\\cite{Horn2004a} with point charges representing ions, and SWM4-NDP~\\cite{Lamoureux2006a} with ions due to Yu {\\em et al.}~\\cite{Yu2010a}, point core+Gaussian vsite (GC+PGV), and point charge + Gaussian vsite and shell (PC+GVS) using ACT.}")
    file.write("\n")
    file.write("\\label{tab:ion_water2}")
    file.write("\n")
    file.write("\\begin{tabular}{lcccccc} \n")
    file.write("\\hline \n")
    file.write("Ion & r$_{min}$ & SAPT & TIP4P-Ew & SWM4-NDP & GC+PGV & PC+GVS\\\\\n")
    file.write("\\hline \n")


    for i in range(len(data['Ion'])):
        file.write(f"{data['Ion'][i]} & {data['r$_{min}$'][i]} & {data['E$_{elec}$ (kJ/mol) SAPT'][i]:.1f} & {data['E$_{elec}$ (kJ/mol) 4W'][i]:.1f} & {data['E$_{elec}$ (kJ/mol) SW'][i]:.1f} & {data['E$_{elec}$ (kJ/mol) ACT$_{GC}$'][i]:.1f}  & {data['E$_{elec}$ (kJ/mol) ACT$_{SS}$'][i]:.1f} \\\\\n")
    file.write("\\hline\n")
    file.write(f"RMSD & & & {rmsd_value1:.1f} & {rmsd_value11:.1f} & {rmsd_value22:.1f} & {rmsd_value2:.1f}\\\\\n")
    file.write(f"MSE & & & {mse_value1:.1f} & {mse_value11:.1f} & {mse_value22:.1f} & {mse_value2:.1f} \\\\\n")

    file.write("\\hline \n")
    file.write("\\end{tabular} \n")
    file.write("\\end{table}")

print("Please check file %s" % file_path)
