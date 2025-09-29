#!/usr/bin/env python3

import os
from run_gaussian import get_mols

debug = False

def extract_columns(mol:str, input_file, skip_lines=7)->list:
    renum = { 
        "ammonium": { "4": 2, "5": 1, "6": 3, "7": 4, "8": 5 },
        "acetate":  { "4": 1, "5": 5, "6": 6, "7": 7, "8": 2, "9": 3, "10": 4 },
        "bromide": { "4": 1 }, "chloride": { "4": 1 }, "fluoride": { "4": 1 },
        "lithium-ion": { "4": 1 }, "sodium-ion": { "4": 1 }, "potassium-ion": { "4": 1 },
        "butanoate": { "4": 3, "5": 2, "6": 4, "7": 1, "8": 6, "9": 7, "10": 5,
                       "11": 8, "12": 9, "13": 10, "14": 11, "15": 12, "16": 13 },
        "ethylammonium": {  "4": 1, "5": 3, "6": 4, "7": 5, "8": 2, "9": 6, "10": 7,
                            "11": 8, "12": 9, "13": 10, "14": 11 },
        "formate": { "4": 1, "5": 3, "6": 4, "7": 2 },
        "guanidine": { "4": 2, "5": 7, "6": 8, "7": 1, "8": 4, "9": 9, "10": 3,
                       "11": 5, "12": 6 },
        "guanidinium": { "4": 2, "5": 7, "6": 8, "7": 1, "8": 3, "9": 5, "10": 6,
                         "11": 4, "12": 9, "13": 10 },
        "imidazole": { "4": 1, "5": 2, "6": 3, "7": 4, "8": 5, "9": 6, "10": 7,
                         "11": 8, "12": 9 },
        "imidazolium": { "4": 1, "5": 6, "6": 2, "7": 7, "8": 3, "9": 8, "10": 4,
                         "11": 9, "12": 5, "13": 10 },
        "1-methylimidazolium": { "4": 3, "5": 8, "6": 2, "7": 7, "8": 1, "9": 13, "10": 5,
                                 "11": 9, "12": 4, "13": 6, "14": 10, "15": 11, "16": 12 },
        "methylammonium": { "4": 1, "5": 3, "6": 4, "7": 5, "8": 2, "9": 6, "10": 7,
                            "11": 8 },
        "propanoate": { "4": 1, "5": 5, "6": 2, "7": 3, "8": 6, "9": 7, "10": 4,
                        "11": 8, "12": 9, "13": 10 },
        "water": { "4": 1, "5": 2, "6": 3 }
    }
    if not mol in renum:
        return None
    if not os.path.exists(input_file):
        return None
    # Create table with right length
    atomq = [-1]*len(renum[mol])
    with open(input_file, 'r') as infile:
        for _ in range(skip_lines):
            next(infile)
        
        for line in infile:
            if line.strip() == 'LOOP':
                break
            
            parts = line.split()
            if len(parts) == 11 and parts[1] != 'DUMM':
                index = parts[0]
                if not index in renum[mol]:
                    print(f"Cannot find index {index} in renum[{mol}]")
                    return None
                atomq[renum[mol][index]-1] = { "name": parts[1], "q": parts[-1], "index": renum[mol][index] }
    return atomq

if __name__ == "__main__":
    mols = get_mols()
    for mol in mols.keys():
        mdir = f"HF_SC/{mol}"
        if os.path.isdir(mdir):
            txml   = f"{mdir}/temp.xml"
            molxml = f"{mdir}/{mol}.xml"
            os.system(f"gauss2molprop -n {mol} -i {mdir}/{mol}.log -o {txml} -basis aug-cc-pvtz")
            
            atomq = {}
            for method in [ "bcc", "resp" ]:
                atomq[method] = extract_columns(mol, f"prepi/{mol}_{method}.prepi")
            with open(molxml, "w") as outf:
                with open(txml) as fd:
                    iatom = 0
                    for line in fd:
                        outf.write(line)
                        if line.find("qMulliken") >= 0:
                            for method in [ "bcc", "resp" ]:
                                if method in atomq and atomq[method]:
                                    m   = method.upper()
                                    if debug:
                                        print(f"mol {mol} iatom {iatom} atomq[method] {atomq[method]}")
                                    myq = atomq[method][iatom]["q"]
                                    outf.write(f"      <q{m}>{myq}</q{m}>\n")
                            iatom += 1
                os.unlink(txml)
            
    os.system("alexandria edit_mp -mp HF_SC/*/*.xml -o hf-aug-cc-pvtz.xml")
