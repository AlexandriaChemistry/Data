#!/usr/bin/env python3

def read_sel(filenm:str)->dict:
    mysel = {}
    with open(filenm, "r") as inf:
        for line in inf:
            words = line.strip().split("|")
            if len(words) == 2:
                mysel[words[0]] = words[1]
    return mysel
    
if __name__ == "__main__":
    prefix = "../../ACTdata/Training/"
    mysel = {}
    for sel in [ 1, 2, 3 ]:
        filenm = prefix + f"ions-nobles-gases{sel}.dat"
        mysel[sel] = read_sel(filenm)
    seltex = "selection.tex"
    with open(seltex, "w") as outf:
        outf.write("\\begin{longtable}{lccc}\n")
        outf.write("\\caption{Data sets used in training. If a cross is set, the particular dimer is part of the training in that set.}\\\\\n")
        outf.write("\\hline\n")
        outf.write("Dimer & Set 1 & Set 2 & Set 3\\\\\n")
        outf.write("\\hline\n")
        for mm in sorted(mysel[1].keys()):
            dimer = mm.replace("#", "-")
            outf.write(f"{dimer}")
            for sel in [ 1, 2, 3 ]:
                if mysel[sel][mm] == "Train":
                    outf.write("& X")
                else:
                    outf.write("& ")
            outf.write("\\\\\n")
        outf.write("\\hline\n")
        outf.write("\\end{longtable}\n")
    print("Please check %s" % seltex)
