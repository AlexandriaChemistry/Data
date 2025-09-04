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
                if len(words) == 3:
                    mydict[words[0]] = { "min": float(words[1]),
                                         "max": float(words[2]) }
    return mydict

def get_nsapt(filenm)->dict:
    nsapt = {}
    with open(filenm, "r") as inf:
        for line in inf:
            if line.startswith("Info: Computing"):
                words = line.strip().split()
                nsapt[words[2]] = words[7]
    return nsapt

if __name__ == "__main__":
    mindist = read_mindist()
    prefix  = "../../ACTdata/Training/"
    mysel   = {}
    Train   = "Train"
    Test    = "Test"
    Total   = "Total"
    # Number in each category
    nnsel   = { Train: [], Test: [], Total: [] }
    allsel = [ 1, 2, 3, 4, 5 ]
    for sel in allsel:
        filenm = prefix + f"ions-nobles-gases{sel}.dat"
        mysel[sel] = read_sel(filenm)
        for tt in [ Train, Test, Total ]:
            nnsel[tt].append(0)
    Nsapt = get_nsapt("logs-MORSEic-Kronecker-EXTRA/Elec-1-C.log")
    seltex = "selection.tex"
    with open(seltex, "w") as outf:
        outf.write("\\begin{longtable}{lccc")
        for c in range(len(allsel)):
            outf.write("c")
        outf.write("}\n")
        outf.write("\\caption{Data sets used in training. If a cross is set, the particular dimer is part of the training in that set. N is the number of data points (SAPT calculations) for this dimer, $r_{min}$ and $r_{max}$ indicate the shortest and longest distances (\\AA) between atoms considered in SAPT calculations and force field training. The final three rows report on the total number of SAPT calculations in each data set.}\\\\\n")
        outf.write("\\hline\n")
        outf.write("Dimer ")
        for sel in allsel:
            outf.write("& Set %d" % sel)
        outf.write(" & N & $r_{min}$ & $r_{max}$ \\\\\n")
        outf.write("\\hline\n")
        for mm in sorted(mysel[1].keys()):
            dimer = mm.replace("#", "-")
            outf.write(f"{dimer}")
            ns = "N/A"
            numsapt = 0
            if mm in Nsapt:
                ns = Nsapt[mm]
                numsapt = int(ns)
            for sel in allsel:
                if mysel[sel][mm] == Train:
                    outf.write("& x")
                    nnsel[Train][sel-1] += numsapt
                    nnsel[Total][sel-1] += numsapt
                else:
                    outf.write("& ")
                    nnsel[Test][sel-1]  += numsapt
                    nnsel[Total][sel-1] += numsapt
            if mm in mindist:
                outf.write(" & %s & %.2f & %.2f\\\\\n" %
                           ( ns, mindist[mm]["min"], mindist[mm]["max"] ) )
            else:
                outf.write(" & &  \\\\\n")
        outf.write("\\hline\n")
        for tt in [ Train, Test, Total ]:
            outf.write("\\# %s " % tt)
            for sel in allsel:
                outf.write(" & %d" % ( nnsel[tt][sel-1] ) )
            outf.write(" & & &\\\\\n")
        outf.write("\\hline\n")
        outf.write("\\end{longtable}\n")
    print("Please check %s" % seltex)
