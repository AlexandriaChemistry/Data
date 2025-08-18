#!/usr/bin/env python3

def read_sel(filenm:str)->dict:
    mysel = {}
    with open(filenm, "r") as inf:
        for line in inf:
            words = line.strip().split("|")
            if len(words) == 2:
                mysel[words[0]] = words[1]
    return mysel

def read_mindist()->dict:
    mydict = {}
    with open("distance.csv", "r") as inf:
        for line in inf:
            if not line.startswith("#"):
                words = line.strip().split(",")
                if len(words) == 2:
                    mydict[words[0]] = float(words[1])
    return mydict

if __name__ == "__main__":
    mindist = read_mindist()
    prefix  = "../../ACTdata/Training/"
    mysel   = {}
    for sel in [ 1, 2, 3 ]:
        filenm = prefix + f"ions-nobles-gases{sel}.dat"
        mysel[sel] = read_sel(filenm)
    seltex = "selection.tex"
    with open(seltex, "w") as outf:
        outf.write("\\begin{longtable}{lcccc}\n")
        outf.write("\\caption{Data sets used in training. If a cross is set, the particular dimer is part of the training in that set. $r_{min}$ indicates the shortest distance (\\AA) between atoms considered in SAPT calculations and force field training.}\\\\\n")
        outf.write("\\hline\n")
        outf.write("Dimer & Set 1 & Set 2 & Set 3 & $r_{min}$ \\\\\n")
        outf.write("\\hline\n")
        for mm in sorted(mysel[1].keys()):
            dimer = mm.replace("#", "-")
            outf.write(f"{dimer}")
            for sel in [ 1, 2, 3 ]:
                if mysel[sel][mm] == "Train":
                    outf.write("& X")
                else:
                    outf.write("& ")
            if mm in mindist:
                outf.write(" & %.1f\\\\\n" % mindist[mm])
            else:
                outf.write(" & \\\\\n")
        outf.write("\\hline\n")
        outf.write("\\end{longtable}\n")
    print("Please check %s" % seltex)
