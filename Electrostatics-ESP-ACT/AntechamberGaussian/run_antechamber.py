#!/usr/bin/env python3

import os, glob
from run_gaussian import get_mols

if __name__ == "__main__":
    mols = get_mols()
    pp   = "prepi"
    os.makedirs(pp, exist_ok=True)
    os.chdir(pp)
    for mol in mols.keys():
        xyz  = f"HF_SC/xyz/{mol}.xyz"
        pdb  = f"HF_SC/{mol}/{mol}.pdb"
        os.system(f"obabel -ixyz {xyz} -opdb -O {pdb}")
        for method in [ "resp", "bcc" ]:
            os.system(f"antechamber -i ../HF_SC/{mol}/{mol}.log -fi gout -rn UNL -o {mol}_{method}.prepi -fo prepi -c {method} -nc {mols[mol]}")
    os.chdir("..")

