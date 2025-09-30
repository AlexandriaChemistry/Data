#!/usr/bin/env python3

import glob, os

os.chdir("HF_SC")
for mol in glob.glob("*"):
    if os.path.isdir(mol) and not mol == "xyz":
        logf = f"{mol}/{mol}.log"
        if os.path.exists(logf):
            xyz = f"../xyz/{mol}.xyz"
            os.system(f"obabel -ig09 {logf} -oxyz -O {xyz}")
