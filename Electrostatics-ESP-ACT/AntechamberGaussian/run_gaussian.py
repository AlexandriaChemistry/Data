#!/usr/bin/env python3

import os

def get_mols()->dict:
    mols = {}
    with open("mols.csv", "r") as inf:
        for line in inf:
            words = line.strip().split(",")
            if len(words) == 2:
                mols[words[0]] = int(words[1])
    return mols

def write_job(mol:str, q:int)->str:
    xyz  = f"../xyz/{mol}.xyz"
    if not os.path.exists(xyz):
        print("No coordinate file %s" % xyz)
        return None
    with open(xyz, "r") as inf:
        lines = inf.readlines()
        com = mol + ".com"
        with open(com, "w") as outf:
            outf.write("""%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            """)
            outf.write("\n%s\n\n" %mol)
            outf.write("%d 1\n" % q)
            for n in range(2, len(lines)):
                outf.write(lines[n])
            outf.write("\n\n")
        job = mol + ".sh"
        with open(job, "w") as outf:
            outf.write("""#!/bin/sh
# Note these commands are cluster-specific
#SBATCH -A naiss2024-3-13
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 72:00:00
            """)
            outf.write("g16 %s\n" % com)
        return job

def run_it():
    mols = get_mols()
    hfdir = "HF_SC"
    os.makedirs(hfdir, exist_ok=True)
    os.chdir(hfdir)
    for mol in mols.keys():
        os.makedirs(mol, exist_ok=True)
        os.chdir(mol)
        # Do not redo a calc if it is already there
        if not os.path.exists(f"{mol}.log"):
            job = write_job(mol, mols[mol])
            if job:
                os.system("sbatch "+job)
        os.chdir("..")
    os.chdir("..")

if __name__ == "__main__":
    run_it()

