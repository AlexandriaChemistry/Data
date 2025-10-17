#!/usr/bin/env python3
import re, os
import json
import subprocess
import sys


log_file_map = {
    "ACM-PG.log": "PC+GS",
    "SWM4.log": "CHARMM Drude"
}

molprops = "../AlexandriaFF/sapt-esp5.xml"
charges = "../AlexandriaFF/hf-aug-cc-pvtz.xml"
selection_file = "../Selection/water+ions.dat"


forcefield_map = {
    "ACM-PG.log": "../AlexandriaFF/PC+GS-elec.xml",
    "SWM4.log": "../ForceFields/CharmmDrude.xml"
}


def extract_rmsd(log_filename):
    rmsd_value = None
    with open(log_filename, "r") as file:
        for line in file:
            if "INDUCTION (kJ" in line:
                parts = line.strip().split()
                try:
                    rmsd_value = round(float(parts[5]), 1)
                except (ValueError, IndexError):
                    print(f"Warning: Could not parse RMSD from line: {line.strip()}")
    if rmsd_value is None:
        print(f"Warning: No RMSD found in {log_filename}")
        rmsd_value = 0.0
    return rmsd_value


for log_filename, ff_xml in forcefield_map.items():
    ff_path = ff_xml
    command = [
        "alexandria", "train_ff", "-nooptimize",
        "-g", log_filename,
        "-sel", selection_file,
        "-mp", molprops,
        "-charges", charges,
        "-ff", ff_path
    ]
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)


    cmd = "grep -v Alexandria INDUCTION.xvg | grep -v Train"
    if "SWM4" in log_filename:
        cmd += " > INDUCTION-SW.xvg"
    else:
        cmd += " > INDUCTION-ACM.xvg"
    os.system(cmd)

rmsd_swm4 = extract_rmsd("SWM4.log")
rmsd_acm = extract_rmsd("ACM-PG.log")


viewxvg_command = [
    "viewxvg", "-f", "INDUCTION-SW.xvg", "INDUCTION-ACM.xvg",
    "-mk", "o", "-res", "-ls", "None", "-color", "green",
    "-alfs", "32", "-tfs", "32", "-lfs", "32",
    "-labels",
    f"CHARMM Drude (RMSD:{rmsd_swm4})", f"PC+GVS (RMSD:{rmsd_acm})",
    "-tickfs", "32", "-color", "crimson", "royalblue",
    "-legend_y", "0.3", "-legend_x", "0.2"
]
print("Running:", " ".join(viewxvg_command))
subprocess.run(viewxvg_command, check=True)
