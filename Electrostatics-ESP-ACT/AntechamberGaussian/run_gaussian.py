#!/usr/bin/env python3

import os

def get_mols()->dict:
    mols = {}
    with open("mols.csv", "r") as inf:
        for line in inf:
            words = line.strip().split(",")
            if len(words) == 3:
                mols[words[0]] = { "q": int(words[1]), "basis": words[2] }
    return mols

def write_job(mol:str, q:int, method:str, basis:str, ncore:int)->str:
    xyz  = f"../../xyz/{mol}.xyz"
    if not os.path.exists(xyz):
        print("No coordinate file %s" % xyz)
        return None
    density = ""
    if method == "MP2":
        density = " density=MP2"
    with open(xyz, "r") as inf:
        lines = inf.readlines()
        com = mol + ".com"
        with open(com, "w") as outf:
            outf.write("""%%mem=12000MB
%%nprocshared=%d
%%mem=%dGb
%%chk=%s.chk
#P %s/aug-cc-pvtz Opt=tight %s Pop=(MK,Hirshfeld) iop(6/33=2) iop(6/42=6)
maxdisk=128GB
            """ % ( ncore, ncore*2, mol, method, density ) )
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
#SBATCH -c %d
#SBATCH -t 72:00:00
            """ % ncore)
            outf.write("g16 %s\n" % com)
        return job

def run_it():
    mols = get_mols()
    for method in [ "HF", "MP2" ]:
        hfdir = f"{method}_SC"
        os.makedirs(hfdir, exist_ok=True)
        os.chdir(hfdir)
        for mol in mols.keys():
            os.makedirs(mol, exist_ok=True)
            os.chdir(mol)
            # Do not redo a calc if it is already there
            if not os.path.exists(f"{mol}.log"):
                job = write_job(mol, mols[mol]["q"], method, mols[mol]["basis"], 16)
                if job:
                    os.system("sbatch "+job)
            os.chdir("..")
        os.chdir("..")

if __name__ == "__main__":
    run_it()

