#!/usr/bin/env python3
import re
import json

log_file_map = {
    "log/walz.log": "eelec-Walz2018a",
    "log/pg.log":   "eelec-ACT-S",
    "log/gc.log":   "eelec-ACT-G",
    "log/jc.log":   "eelec-CT"
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
    "potassium#bromide":  {"rmin": 2.821, "eelec-SAPT":  -538.38}
}

def extract_log_value(log_file, molecule, qm_value):
    try:
        with open(log_file, 'r') as file:
            content = file.read()
        match = re.search(fr"Name: {re.escape(molecule)}.*?(?:\n\s*\n|$)", content, re.DOTALL)
        if not match:
            return None

        section = match.group(0)
        for line in section.strip().splitlines():
            if str(qm_value) in line:
                tokens = line.split()
                try:
                    idx = tokens.index(str(qm_value))
                    return float(tokens[idx + 1]) if idx + 1 < len(tokens) else None
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
ø
