#!/usr/bin/env python3
from tabulate import tabulate
import numpy as np
import json



with open('sapt2_data.json', 'r') as json_file:
    sapt2 = json.load(json_file)


def rmsd(data, ref_key1, ref_key2):
    ref1_values = []
    ref2_values = []
    for ion_pair, values in data.items():
        if ref_key1 in values and ref_key2 in values:
            ref1_values.append(values[ref_key1])
            ref2_values.append(values[ref_key2])
    return np.sqrt(np.mean((np.array(ref1_values) - np.array(ref2_values))**2))


def mse(data, ref_key1, ref_key2):
        ref1_values = []
        ref2_values = []
        for ion_pair, values in data.items():
            if ref_key1 in values and ref_key2 in values:
                ref1_values.append(values[ref_key1])
                ref2_values.append(values[ref_key2])
        return np.mean(np.subtract(ref1_values,ref2_values))


def dict_to_list_of_lists(data):
    rows = []
    for key, value in data.items():
        rows.append([key.replace("-",""), value.get("rmin", "N/A"),  f"{value.get('eelec-SAPT', 'N/A'):.1f}", f"{value.get('eelec-CT', 'N/A'):.1f}", f"{value.get('eelec-Walz2018a', 'N/A'):.1f}",  f"{value.get('eelec-ACT-G', 'N/A'):.1f}", f"{value.get('eelec-ACT-S', 'N/A'):.1f}"])
    return rows


headers = ["Ions", "r$_{min}$", "SAPT", "PC", "Walz {\\em et al.}", "GC+PGV", "PC+GVS"]


table_string= tabulate(dict_to_list_of_lists(sapt2), headers, tablefmt="latex_raw")


caption= "\\caption{Minimum energy distance (\\AA) between ions and electrostatic energies from SAPT2+(CCD)$\\delta$MP2 with the aug-cc-pVTZ basis set, for point charges (PC), the Walz {\\em et al.} model with a Gaussian charge distribution~\\cite{Walz2018a}, and the ACT models GC+PVG and PC+GVS consisting of a point charge and virtual site with a Gaussian charge (K$^+$ and halide ions only) and a Drude particle with a Gaussian charge (PC+GVS only). The RMSD and MSE were calculated with respect to the SAPT2+(CCD)$\\delta$MP2 with the aug-cc-pVTZ basis set electrostatic energy. Note that the level of theory is different from Tables~\\ref{tab:sapt100} and~\\ref{tab:sapt10} and the results cannot be compared directly. In addition, the results in Tables~\\ref{tab:sapt100} and~\\ref{tab:sapt10} are from fitting models to the ESP, whereas in this table training was done on SAPT data as indicated above.}"


label = "\\label{tab:sapt2_ions}"



rmsd_value1 = rmsd(sapt2, "eelec-CT", "eelec-SAPT")
rmsd_value2 = rmsd(sapt2, "eelec-ACT-S", "eelec-SAPT")
rmsd_value3 = rmsd(sapt2, "eelec-ACT-G", "eelec-SAPT")
rmsd_value6 = rmsd(sapt2, "eelec-Walz2018a", "eelec-SAPT")


mse_value1 = mse(sapt2, "eelec-CT", "eelec-SAPT")
mse_value2 = mse(sapt2, "eelec-ACT-S", "eelec-SAPT")
mse_value3 = mse(sapt2, "eelec-ACT-G", "eelec-SAPT")
mse_value6 = mse(sapt2, "eelec-Walz2018a", "eelec-SAPT")


#table = '\n'.join(line for line in table_string.split('\n')[3:-2])



table_with_rmsd = table_string + f"\nRMSD & & & {rmsd_value1:.1f} & {rmsd_value6:.1f} & {rmsd_value3:.1f} & {rmsd_value2:.1f} \\\\"
table_with_rmsd = table_with_rmsd + f"\nMSE & & & {mse_value1:.1f} & {mse_value6:.1f} & {mse_value3:.1f} & {mse_value2:.1f} \\\\"

file_path = "Ions-sapt2-JC-Walz2018a-ACT.tex"


with open(file_path, "w") as file:
    file.write("\\begin{table}[ht]\n")
    file.write("\\centering\n")
    file.write(caption + "\n" + label + "\n")
    for line in table_with_rmsd:
        file.write(line)
    file.write("\n")
    file.write("\\end{tabular}\n")
    file.write("\\end{table}\n")
    
print("Please fix %s manually to correct output." % file_path)




