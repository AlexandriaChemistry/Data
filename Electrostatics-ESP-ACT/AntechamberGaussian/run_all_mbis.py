#!/usr/bin/env python3

import os
from run_gaussian import get_mols

def run_it(compound:str, charge:int, basis:str):
    os.chdir(compound)
    method = "CCSD"
    if not (os.path.exists("scf.json") and os.path.exists(f"{method}.json")):
        os.system(f"sbatch ../run_mbis.py {compound} {charge} {method} {basis}")
    os.chdir("..")

if __name__ == "__main__":
    mols = get_mols()
    os.chdir("MBIS")
    for mol in mols:
        os.makedirs(mol, exist_ok=True)
        run_it(mol, mols[mol]["q"], mols[mol]["basis"])
