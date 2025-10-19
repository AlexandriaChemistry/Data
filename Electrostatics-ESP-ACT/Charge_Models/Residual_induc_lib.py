#!/usr/bin/env python3
import os
import subprocess
import matplotlib.pyplot as plt

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
    command = [
        "alexandria", "train_ff", "-nooptimize",
        "-g", log_filename,
        "-sel", selection_file,
        "-mp", molprops,
        "-charges", charges,
        "-ff", ff_xml
    ]
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)

    output_xvg = "INDUCTION-SW.xvg" if "SWM4" in log_filename else "INDUCTION-ACM.xvg"
    os.system(f"grep -v Alexandria INDUCTION.xvg | grep -v Train > {output_xvg}")

rmsd_swm4 = extract_rmsd("SWM4.log")
rmsd_acm = extract_rmsd("ACM-PG.log")

def read_xvg(filename):
    x, y = [], []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@", "&")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
            except ValueError:
                print(f"Warning: Could not parse line in {filename}: {line}")
    return x, y

x_swm4, y_swm4 = read_xvg("INDUCTION-SW.xvg")
x_acm, y_acm = read_xvg("INDUCTION-ACM.xvg")

plt.figure(figsize=(10, 6))

plt.scatter(
    x_swm4, y_swm4,
    color='royalblue', alpha=0.7,
    s=120, edgecolor='black', linewidth=0.8,
    marker='o',
    label=f"CHARMM Drude"
)

plt.scatter(
    x_acm, y_acm,
    color='crimson', alpha=0.7,
    s=120, edgecolor='black', linewidth=0.8,
    marker='s',
    label=f"PC+GS"
)

plt.xlabel("Induction (kJ/mol)", fontsize=32)
plt.ylabel("Residual (kJ/mol)", fontsize=32)

plt.xticks(fontsize=32)
plt.yticks(fontsize=32)

#plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.4)

plt.legend(
    loc='lower right', fontsize=24,
    frameon=True, shadow=True, facecolor='white'
)

plt.text(0.02, 0.95, f"SWM4 RMSD: {rmsd_swm4}", transform=plt.gca().transAxes,
         fontsize=20, color='royalblue', weight='bold')
plt.text(0.02, 0.85, f"ACM RMSD: {rmsd_acm}", transform=plt.gca().transAxes,
         fontsize=20, color='crimson', weight='bold')

plt.tight_layout()
plt.savefig("induction-water-ions.pdf", dpi=300)
plt.show()

