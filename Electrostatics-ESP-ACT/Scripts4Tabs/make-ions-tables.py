#!/usr/bin/env python

import json, math, os, sys

debug = True

def logToCoulomb(logfn:str)->float:
    with open(logfn, "r") as inf:
        for line in inf:
            if "Info: Interaction energy COULOMB:" in line:
                words = line.split()
                if len(words) != 5:
                    sys.exit("Expected 5 words on line '%s'" % line.strip())
                try:
                    return float(words[4])
                except ValueError:
                    sys.exit("Cannot interpret line '%s'" % line.strip())

def logToInduction(logfn:str)->float:
    einduc = 0
    with open(logfn, "r") as inf:
        for line in inf:
            if ( "Info: Interaction energy INDUCTION:" in line or
                 "Info: Interaction energy INDUCTIONCORRECTION:" in line):
                words = line.split()
                if len(words) != 5:
                    sys.exit("Expected 5 words on line '%s'" % line.strip())
                try:
                    einduc += float(words[4])
                except ValueError:
                    sys.exit("Cannot interpret line '%s'" % line.strip())
    return einduc

def xvgToEner(xvgfn:str)->float:
    with open(xvgfn, "r") as inf:
        myset = -1
        for line in inf:
            if "COULOMB" in line and "legend" in line:
                words = line.split()
                myset = int(words[1][1:])
            elif not ("#" in line) and not ("@" in line):
                words = line.split()
                try:
                    if len(words) > myset+1:
                        return float(words[myset+1])
                    else:
                        continue
                        #sys.exit("Stange line '%s'" % line.strip())
                except ValueError:
                    sys.exit("Cannot interpret line '%s'" % line.strip())

def add_calcs(mydata:list, models:dict):
    tdir = "simulation_output"
    os.makedirs(tdir, exist_ok=True)
    for m in models.keys():
        for dim in range(len(mydata)):
            mylog = tdir + "/" + mydata[dim]["name"] + "-" + m + ".log"
            mytrj = mylog[:-3] + "pdb"
            myxvg = mylog[:-3] + "xvg"
            myconf = ( "Conformations/%s.sdf" % mydata[dim]["name"])
            if not os.path.exists(myconf):
                print("File %s is missing" % myconf)
                continue
            mycmd = ("alexandria simulate -ff %s -charges %s  -f %s  -sp -g %s -o %s -e %s" %
                     ( models[m]["ff"], models[m]["mp"], myconf, mylog, mytrj, myxvg ) )
            if "qtype" in models[m]:
                mycmd += ( " -qqm %s" % models[m]["qtype"] )
            print(mycmd)
            os.system(mycmd)
            mydata[dim][m] = logToCoulomb(mylog)

            if not debug:
                for myfn in [ mylog, mytrj,myxvg ]:
                    if os.path.exists(myfn):
                        os.unlink(myfn)

def wtable(file_path:str, models:dict, mydata:list, caption:str, label:str):
    with open(file_path, "w") as file:
        file.write("\\begin{table}[ht]\n")
        file.write("\\centering\n")
        file.write("\\caption{%s}\n" % caption)
        file.write("\\label{%s}\n" % label)
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
                if not m in mydata[i] or not mydata[i][m]:
                    file.write("& - ")
                    diff = 0
                else:
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

def water_ions():
    mymp  = "../AlexandriaFF/hf-aug-cc-pvtz.xml"
    models = { "TIP4P-Ew": { "ff": "../ForceFields/TIP4PEW-JC.xml", "mp": mymp },
               "CHARMM": { "ff": "../ForceFields/CharmmDrude.xml", "mp": mymp },
               "GC+PGV": { "ff": "../AlexandriaFF/coul-gv.xml", "mp": mymp },
               "PC+GVS": { "ff": "../AlexandriaFF/all-pg.xml", "mp": mymp } }
    newfn = "data-water-ions.json"
    if not os.path.exists(newfn):
        sys.exit("Cannot find %s" % newfn)
    with open(newfn, "r") as inf:
        mydata = json.load(inf)

    add_calcs(mydata, models)

    file_path = "ion-water-SAPT2-TIP4Pew-ACT4S.tex"
    caption = "\\textbf{Water-ion electrostatic energies at distances close to their energy minimum.} Minimum energy distance (\\AA) between ions and water oxygen/hydrogen from Experiment (ref.~\\citenum{Heyrovska2006a}), and minimized water dimer (ref.~\\citenum{temelso2011benchmark}). Electrostatic energies are reported in kJ/mol from the SAPT2+(CCD)-$\\delta$MP2 method with an aug-cc-pVTZ basis set, TIP4P-Ew~\\cite{Horn2004a} with point charges representing ions, for the CHARMM drude model of water (SWM4-NDP~\\cite{Lamoureux2006a}) with ions due to Yu {\\em et al.}~\\cite{Yu2010a}, as well as point core+Gaussian vsite (GC+PGV), and point charge + Gaussian vsite and shell (PC+GVS) derived here using ACT."
    label = "tab:ion_water2"
    wtable(file_path, models, mydata, caption, label)
    
def water_ions_induction():
    mymp  = "../AlexandriaFF/hf-aug-cc-pvtz.xml"
    models = { "CHARMM": { "ff": "../ForceFields/CharmmDrude.xml", "mp": mymp },
               "PC+GVS": { "ff": "../AlexandriaFF/all-pg.xml", "mp": mymp } }
    newfn = "data-water-ions-induction.json"
    if not os.path.exists(newfn):
        sys.exit("Cannot find %s" % newfn)
    with open(newfn, "r") as inf:
        mydata = json.load(inf)

    tdir = "simulation_output"
    for i in range(len(mydata)):
        for m in models.keys():
            mylog = tdir + "/" + mydata[i]["name"] + "-" + m + ".log"
            if os.path.exists(mylog):
                mydata[i][m] = logToInduction(mylog)
            else:
                mydata[i][m] = None

    file_path = "ion-water-induction.tex"
    caption = "\\textbf{Water-ion induction energies at distances close to their energy minimum.} Minimum energy distance (\\AA) between ions and water oxygen/hydrogen from Experiment (ref.~\\citenum{Heyrovska2006a}), and minimized water dimer (ref.~\\citenum{temelso2011benchmark}). Induction energies are reported in kJ/mol from the SAPT2+(CCD)-$\\delta$MP2 method with an aug-cc-pVTZ basis set, for the CHARMM drude model of water (SWM4-NDP~\\cite{Lamoureux2006a}) with ions due to Yu {\\em et al.}~\\cite{Yu2010a}, as well the point charge + Gaussian vsite and shell (PC+GVS) derived here using ACT."
    label = "tab:ion_inducs"
    wtable(file_path, models, mydata, caption, label)
    
def ac_mt_gaff():
    mymp  = "../AlexandriaFF/hf-aug-cc-pvtz.xml"
    epg   = "../AlexandriaFF/esp-paper-gaussian.xml"
    models = { "RESP": { "ff": "../AlexandriaFF/coul-p.xml", "mp": epg, "qtype": "qRESP" },
               "BCC": { "ff": "../AlexandriaFF/coul-p.xml", "mp": epg, "qtype": "qBCC" },
               "GC+PGV": { "ff": "../AlexandriaFF/coul-gv.xml", "mp": mymp },
               "PC+GVS": { "ff": "../AlexandriaFF/all-pg.xml", "mp": mymp } }
    newfn = "data-sc-ions.json"
    if not os.path.exists(newfn):
        sys.exit("Cannot find %s" % newfn)
    with open(newfn, "r") as inf:
        mydata = json.load(inf)
    add_calcs(mydata, models)

    file_path = "AC-MA-IONS-GAFF.tex"
    caption = "Electrostatic energy (kJ/mol) between alkali ions, halides or water (oxygen) and amino acid side chain analogs, formate (oxygen), acetate (oxygen), methylammonium (nitrogen), ethylammonium (nitrogen) from SAPT2+(CCD)$\\delta$MP2/aug-cc-pVTZ, and charges determined using either RESP~\\cite{Bayly1993a} or BCC~\\cite{Jakalian2000a}, as well as two models generated using the ACT."
    label = "tab:ac_ma_ions"
    wtable(file_path, models, mydata, caption, label)
    
if __name__ == "__main__":
    water_ions()
    water_ions_induction()
#    ac_mt_gaff()
