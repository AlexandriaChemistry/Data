#!/usr/bin/env python3
import re
import json
import subprocess

log_file_map = {
    "Walz.log": "Walz2018a",
    "ACM-all-PG.log":   "PC+GVS",
    "ACM-elec-GV.log":   "GC+PGV",
    "TIP4PEW-JC.log":   "JC",
    "SWM4.log": "SWM4",
    "RESP.log":   "RESP",
    "BCC.log":   "BCC"

}

molecule_data = {
    "lithium#fluoride":   {"rmin": 1.564, "eelec-SAPT":  -850.53},
    "lithium#chloride":   {"rmin": 2.021, "eelec-SAPT":  -656.85},
    "lithium#bromide":    {"rmin": 2.170, "eelec-SAPT":  -609.54},
    "sodium#fluoride":    {"rmin": 1.926, "eelec-SAPT":  -750.62},
    "sodium#chloride":    {"rmin": 2.361, "eelec-SAPT":  -607.90},
    "sodium#bromide":     {"rmin": 2.502, "eelec-SAPT":  -572.68},
    "potassium#fluoride": {"rmin": 2.171, "eelec-SAPT":  -719.80},
    "potassium#chloride": {"rmin": 2.667, "eelec-SAPT":  -570.19},
    "potassium#bromide":  {"rmin": 2.821, "eelec-SAPT":  -538.38},
    "water#lithium":      {"rmin": 2.04,  "eelec-SAPT":  -111.55},
    "water#sodium":       {"rmin": 2.39,  "eelec-SAPT":   -88.71},
    "water#potassium":    {"rmin": 2.80,  "eelec-SAPT":   -68.50},
    "water#fluoride":     {"rmin": 2.64,  "eelec-SAPT":  -146.66},
    "water#chloride":     {"rmin": 3.17,  "eelec-SAPT":   -76.85},
    "water#bromide":      {"rmin": 3.35,  "eelec-SAPT":   -60.49},
    "water#water":        {"rmin": 2.917, "eelec-SAPT":   -31.59},
    "formate#lithium":    {"rmin": 1.848, "eelec-SAPT":  -734.50},
    "formate#sodium":     {"rmin": 2.175, "eelec-SAPT":  -673.98},
    "formate#potassium":  {"rmin": 2.521, "eelec-SAPT":  -604.73},
    "formate#water":      {"rmin": 3.109, "eelec-SAPT":  -136.61},
    "acetate#lithium":    {"rmin": 2.032, "eelec-SAPT":  -624.68},
    "acetate#sodium":     {"rmin": 2.031, "eelec-SAPT":  -646.24},
    "acetate#potassium":  {"rmin": 2.519, "eelec-SAPT":  -574.45},
    "acetate#water":      {"rmin": 3.098, "eelec-SAPT":  -142.25},
    "methylammonium#fluoride": {"rmin": 2.711, "eelec-SAPT": -458.91},
    "methylammonium#chloride": {"rmin": 2.986, "eelec-SAPT": -452.92},
    "methylammonium#bromide":  {"rmin": 3.331, "eelec-SAPT": -434.94},
    "methylammonium#water":    {"rmin": 2.65,  "eelec-SAPT": -100.75},
    "ethylammonium#fluoride":  {"rmin": 2.794, "eelec-SAPT": -440.31},
    "ethylammonium#chloride":  {"rmin": 3.259, "eelec-SAPT": -399.66},
    "ethylammonium#bromide":   {"rmin": 3.391, "eelec-SAPT": -385.98},
    "ethylammonium#water":     {"rmin": 2.651, "eelec-SAPT":  -99.48}
}


molprops = "../AlexandriaFF/sapt-esp.xml"
selection_file = "../Selection/ac-total.dat"
forcefield_map = {
    "Walz.log": "Walz2018a.xml",
    "TIP4PEW-JC.log": "TIP4PEW-JC.xml",
    "SWM4.log": "CharmmDrude.xml"
}

for log_filename, ff_xml in forcefield_map.items():
    ff_path = f"../ForceFields/{ff_xml}"
    command = [
        "alexandria", "train_ff", "-nooptimize",
        "-g", log_filename,
        "-sel", selection_file,
        "-mp", molprops,
        "-ff", ff_path
    ]
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)

def extract_log_value(log_file, molecule, qm_value):
    try:
        with open(log_file, 'r') as file:
            content = file.read()

        match = re.search(fr"Name: {re.escape(molecule)}.*", content, re.DOTALL)
        if not match:
            return None

        section = match.group(0)
        for line in section.strip().splitlines():
            tokens = line.split()
            for i, token in enumerate(tokens):
                try:
                    val = float(token)
                    if abs(val - qm_value) < 0.01 and i + 1 < len(tokens):
                        return float(tokens[i + 1])
                except ValueError:
                    continue
    except Exception as e:
        print(f"Error reading {log_file}: {e}")
    return None


results = {}
for mol, props in molecule_data.items():
    entry = props.copy()
    qm = props["eelec-SAPT"]
    entry["QM"] = qm

    for logfile, label in log_file_map.items():
        value = extract_log_value(logfile, mol, qm)
        entry[label] = value

    results[mol] = entry

with open('sapt2_data.json', 'w') as f:
    json.dump(results, f, indent=4)

print("Data saved to sapt2_data.json")
