#!/usr/bin/env python3

import os, sys

def read_log(filenm:str)->dict:
    mylist = { "Test": {}, "Train": {} }
    with open(filenm, "r") as inf:
        for line in inf:
            for ener in [ "COULOMB", "INDUCTION", "DISPERSION", "EXCHANGE", "EPOT" ]:
                if line.startswith(f"{ener} (kJ/mol)"):
                    words = line.strip().split()
                    rmsd = float(words[5])
                    dset = words[8].split("-")[1]
                    mylist[dset][ener] = { "RMSD": rmsd, "N": int(words[2]) }

    return mylist

if __name__ == "__main__":
    viewflags = "-lfs 32 -alfs 32 -tickfs 28"
    prefix = "../../ACTdata/Training/ions-nobles-gases/TT2b/allhh/pol/vsite/"
    mylist = {}
    allfns = [ "800elec3", "800elec4", "800elec5" ]
    rename = { "COULOMB": "Electrostatics" }
    for fn in allfns:
        mylist[fn] = {}
        for ff in [ "Elec", "Inter" ]:
            logfn = prefix + fn + f"/train-{ff}.log"
            if os.path.exists(logfn):
                allset = read_log(logfn)
                for dset in allset:
                    if not dset in mylist[fn]:
                        mylist[fn][dset] = {}
                    for ener in allset[dset]:
                        newener = ener
                        if ener in rename:
                            newener = rename[ener]
                        mylist[fn][dset][newener] = allset[dset][ener]
    texfn = "stats.tex"
    print(mylist)
    with open(texfn, "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{Root mean square deviation from SAPT energies (kJ/mol) for three different selections (Table S1) separated in Train and Test sets.}\n")
        outf.write("\\begin{tabular}{lcccc}\n")
        outf.write("\\hline\n")
        outf.write("Term & & Set 1 & Set 2 & Set 3\\\\\n")
        outf.write("\\hline\n")
        for dset in [ "Train", "Test" ]:
            for ener in mylist["800elec3"][dset]:
                if "Train" == dset:
                    outf.write("%s & %s" % ( ener, dset ))
                else:
                    outf.write(" & %s" % ( dset ))
                for fn in allfns:
                    outf.write(" & %.1f" % ( mylist[fn][dset][ener]["RMSD"]) )
                outf.write("\\\\\n")
        outf.write("\\hline\n")
        outf.write("\\end{tabular}\n")
        outf.write("\\end{table}\n")
    print("Please check %s" % texfn)
