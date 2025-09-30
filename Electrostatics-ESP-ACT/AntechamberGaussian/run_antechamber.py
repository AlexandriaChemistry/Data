#!/usr/bin/env python3

import os, glob
from run_gaussian import get_mols

if __name__ == "__main__":
    mols = get_mols()
    pp   = "prepi"
    os.makedirs(pp, exist_ok=True)
    os.chdir(pp)
    for qm in [ "HF", "MP2" ]:
        for mol in mols.keys():
            for method in [ "resp", "bcc" ]:
                os.system(f"antechamber -i ../{qm}_SC/{mol}/{mol}.log -fi gout -o {mol}_{qm}_{method}.prepi -fo prepi -c {method} -nc {mols[mol]['q']} -j 4 -pf yes")
    os.chdir("..")

