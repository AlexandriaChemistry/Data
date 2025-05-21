#!/usr/bin/env python

import json, math, os, sys

debug = False

models = { "TIP4P-Ew": "../ForceFields/TIP4PEW-JC.xml",
           "SWM4-NDP": "../ForceFields/SWM4-NDP-ions.xml",
           "GC+PGV": "../AlexandriaFF/coul-gv.xml",
           "PC+GVS": "../AlexandriaFF/all-pg.xml" }

def get_data()->list:
    mydata = []
    newfn = "data-water-ions.json"
    if not os.path.exists(newfn):
        sys.exit("Cannot find %s" % newfn)
    with open(newfn, "r") as inf:
        mydata = json.load(inf)

    return mydata

def add_calcs(mydata:list):
    mymp  = "../AlexandriaFF/hf-aug-cc-pvtz.xml"
    for m in models.keys():
        for dim in range(len(mydata)):
            mylog = mydata[dim]["name"] + "-" + m + ".log"
            os.system("alexandria simulate -ff %s -charges %s  -f Conformations/%s.sdf -minimize -g %s" %
                      ( models[m], mymp, mydata[dim]["name"], mylog ) )
            with open(mylog, "r") as inf:
                for line in inf:
                    coul = 'COULOMB'
                    if coul in line and not "Force" in line and not "Interaction" in line:
                        words = line.split()
                        try:
                            value = float(words[2][:-1])
                            mydata[dim][m] = value
                            break
                        except ValueError:
                            sys.exit("Cannot interpret line '%s'" % line.strip())
            if not debug:
                os.unlink(mylog)

if __name__ == "__main__":
    mydata = get_data()
    add_calcs(mydata)
    if debug:
        print(mydata)
    file_path = "ion-water-SAPT2-TIP4Pew-ACT4S.tex"
    with open(file_path, "w") as file:
        file.write("\\begin{table}[ht]\n")
        file.write("\\centering\n")
        file.write("\\caption{\\textbf{Water-ion energies at their energy minimum.} Minimum energy distance (\\AA) between ions and water oxygen/hydrogen from Experiment (ref.~\\citenum{Heyrovska2006a}), and minimized water dimer (ref.~\\citenum{temelso2011benchmark}). Electrostatic energies are reported in kJ/mol from the SAPT2+(CCD)-$\\delta$MP2 method with an aug-cc-pVTZ basis set, TIP4P-Ew~\\cite{Horn2004a} with point charges representing ions, and SWM4-NDP~\\cite{Lamoureux2006a} with ions due to Yu {\\em et al.}~\\cite{Yu2010a}, point core+Gaussian vsite (GC+PGV), and point charge + Gaussian vsite and shell (PC+GVS) using ACT.}")
        file.write("\n")
        file.write("\\label{tab:ion_water2}")
        file.write("\n")
        file.write("\\begin{tabular}{lcccccc} \n")
        file.write("\\hline \n")
        file.write("Ion & r$_{min}$ & SAPT ")
        for m in models.keys():
            file.write(" & %s " % m)
        file.write("\\\\\n")
        file.write("\\hline \n")

        rmsd = {}
        mse  = {}
        for m in models.keys():
            rmsd[m] = 0
            mse[m]  = 0
        for i in range(len(mydata)):
            file.write("%s & %g & %.1f " %
                       ( mydata[i]["latex"], mydata[i]["rmin"], mydata[i]["sapt2"] ) )
            for m in models.keys():
                file.write(" & %.1f" % mydata[i][m])
                diff = mydata[i][m] - mydata[i]["sapt2"]
                rmsd[m] += diff**2
                mse[m]  += diff
            file.write("\\\\\n")
            
        file.write("\\hline\n")
        file.write("RMSD & &")
        for m in models.keys():
            file.write(" & %.1f " % math.sqrt(rmsd[m]/len(mydata)) )
        file.write("\\\\\n")
        file.write("MSE & &")
        for m in models.keys():
            file.write(" & %.1f " % ( mse[m]/len(mydata) ) )
        file.write("\\\\\n")
        file.write("\\hline \n")
        file.write("\\end{tabular} \n")
        file.write("\\end{table}")

    print("Please check file %s" % file_path)
