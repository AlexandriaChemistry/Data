#!/usr/bin/env python3

import os, glob

if __name__ == "__main__":
    mc   = "HF_SC/mols.csv"
    mols = {}
    with open(mc, "r") as inf:
        for line in inf:
            words = line.strip().split(",")
            if len(words) == 2:
                mols[words[0]] = int(words[1])

    pp = "prepi"
    os.makedirs(pp, exist_ok=True)
    os.chdir(pp)
    for mol in mols.keys():
        for method in [ "resp", "bcc" ]:
            os.system(f"antechamber -i ../HF_SC/{mol}/{mol}.log -fi gout -o {mol}_{method}.prepi -fo prepi -c {method} -nc {mols[mol]}")
    os.chdir("..")

