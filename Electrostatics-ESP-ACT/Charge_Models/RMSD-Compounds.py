#!/usr/bin/env python3

import os, sys, math
import numpy as np
import holoviews as hv
from holoviews import opts
hv.extension('matplotlib')

monomers = [ "water", "lithium", "potassium", "sodium",
             "fluoride", "chloride", "bromide",
             "formate", "acetate", "propanoate", "butanoate",
             "ammonium", "methylammonium", "ethylammonium",
             "guanidinium", "imidazolium" ]

def get_train()->list:
    mylist = []
    with open("../Selection/ac-train.dat", "r") as inf:
        for line in inf:
            ww = line.strip().split("|")
            mylist.append(ww[0])
    return mylist
        
def run_one(filenm: str) -> dict:
    tempf = "rmsd.txt"
    os.system("grep 'COULOMB RMSD' %s > %s" % (filenm, tempf))
    
    mydict = {}
    with open(tempf, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if len(words) == 14:
                mydict[words[13]] = {"RMSD": words[2], "N": words[12]}
    return mydict


def run_legacy(filenm: str) -> dict:
    tempf = "rmsd_legacy.txt"
    os.system("grep 'COULOMB RMSD' %s > %s" % (filenm, tempf))
    
    mydict = {}
    with open(tempf, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if len(words) == 14:
                mydict[words[13]] = words[2]
    return mydict

def make_plot(model:dict):
    newdata = []
    for i in monomers:
        for j in monomers:
            # Making strings out of the equation to keep the
            # boxes equal in size.
            ij = f"{i}#{j}"
            if not ij in model:
                ij = f"{j}#{i}"
            rmsd = None
            if ij in model:
                rmsd = float(model[ij]["RMSD"])
            newentry = ( i, j, rmsd )
            newdata.append(newentry)
            if i != j:
                newentry = ( j, i, rmsd )
                newdata.append(newentry)
    xx = hv.Dimension("xx", values=monomers)
    yy = hv.Dimension("yy", values=monomers)
    return hv.HeatMap(newdata).sort().aggregate(function=np.min)

if __name__ == "__main__":
    alldata = {}
    mydirs = {
        "ESP": "ESP_MP2.log",
        "BCC": "BCC_MP2.log",
        "RESP": "RESP_MP2.log",
        #"MBIS": "MBIS_MP2.log",
        "MBIS-S": "MBIS-S_MP2.log",
        "PC+GV4x": "PC+GV-esp4_MP2.log",
        "PC+SV4x": "PC+SV-esp4_MP2.log",
        "PC": "PC-elec_MP2.log",
        "GC": "GC-elec_MP2.log",
        "SC": "SC-elec_MP2.log",
        "PC+GV4": "PC+GV-elec_MP2.log",
        "PC+SV4": "PC+SV-elec_MP2.log",
        "PC+GS4": "PC+GS-elec_MP2.log",
        #"Mulliken": "Mulliken_MP2.log",
        #"Hirshfeld": "Hirshfeld_MP2.log",
        #"ESP": "ESP_MP2.log",
        #"CM5": "CM5_MP2.log",
    }
    for mydir in mydirs.keys():
        alldata[mydir] = run_one(mydirs[mydir])
        hm = make_plot(alldata[mydir]).opts(show_values=False, title=mydir, xlabel="", xrotation=90,  ylabel="", clabel="RMSD (kJ/mol)", colorbar=True, fontsize='x-large', clim=(0, 40))
        pdfname = ( "heatmap-%s" % ( mydir ))
        hv.save(hm, filename=pdfname, fmt='pdf')

    train = get_train()
    texfn = "rmsdtable.tex"
    with open(texfn, "w") as outf:
        outf.write("\\begin{longtable}{lccccccccccccc}\n")
        outf.write("\\caption{Root mean square deviation (kJ/mol) from SAPT2+(CCD)$\\delta$MP2 electrostatics per compound dimer for some of the ACT models and widely-used models. N is the number of conformations of each dimer used. Compound dimers used in training are printed in {\\bf bold font}.}\\\\\n")
        outf.write("\\label{tab:saptrmsd}\\\\\n")
        outf.write("\\hline\n")
        outf.write("Dimer & N ")
        for md in alldata.keys():
            outf.write(" & \\rotatebox{90}{%s}" % md)
        outf.write("\\\\\n")
        outf.write("\\hline\n")
        for dimer in sorted(alldata["PC"].keys()):
            ddd = dimer.replace("#", "-")
            if dimer in train:
                ddd = ("{\\bf %s}" % ddd)
            outf.write("%s & %s " % (ddd, alldata["PC"][dimer]["N"]))
            for md in alldata.keys():
                if dimer in alldata[md]:
                    outf.write(" & %.1f" % float(alldata[md][dimer]["RMSD"]))
                else:
                    outf.write(" & ")
            outf.write("\\\\\n")
        outf.write("\\hline\n")
        outf.write("\\end{longtable}\n")
    print("Please check %s" % texfn)

