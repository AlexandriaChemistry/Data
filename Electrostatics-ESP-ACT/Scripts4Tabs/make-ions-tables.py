#!/usr/bin/env python

import json, math, os, sys

debug = True
mymp  = "../AlexandriaFF/MP2-aug-cc-pvtz.xml"

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
            mycmd = ("alexandria simulate -ff %s  -f %s  -sp -g %s -o %s -e %s" %
                     ( models[m]["ff"], myconf, mylog, mytrj, myxvg ) )
            if "mp" in models[m] and not models[m]["mp"] is None:
                mycmd += ( " -charges %s " % models[m]["mp"])
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
        file.write("\\begin{tabular}{lc")
        for c in range(len(mydata)):
            file.write("c")
        file.write("} \n")
        file.write("\\hline \n")
        file.write("Ion & r & SAPT ")
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
            myname = mydata[i]["name"].replace("#", "-")
            file.write("%s & %g & %.1f " %
                       ( myname, mydata[i]["rmin"], mydata[i]["sapt2"] ) )
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


def water_ions():
    models = { "TIP4P-Ew": { "ff": "../ForceFields/TIP4PEW-JC.xml", "mp": mymp },
               "MBIS-S": { "ff": "../AlexandriaFF/P+S_updated.xml", "mp": None },
               "CHARMM": { "ff": "../ForceFields/CharmmDrude.xml", "mp": mymp },
               "PC+GV": { "ff": "../AlexandriaFF/PC+GV-elec.xml", "mp": mymp },
               "PC+GS": { "ff": "../AlexandriaFF/PC+GS-elec.xml", "mp": mymp } }
    newfn = "data-water-ions.json"
    if not os.path.exists(newfn):
        sys.exit("Cannot find %s" % newfn)
    with open(newfn, "r") as inf:
        mydata = json.load(inf)

    add_calcs(mydata, models)

    file_path = "ion-water-SAPT2-TIP4Pew-ACT4S.tex"
    caption = "Water-ion electrostatic energies at distances close to their energy minimum. Distance r (\\AA) between ions and water oxygen/hydrogen from Experiment (ref.~\\citenum{Heyrovska2006a}), and minimized water dimer (ref.~\\citenum{temelso2011benchmark}). Electrostatic energies are reported in kJ/mol from the SAPT2+(CCD)-$\\delta$MP2 method with an aug-cc-pVTZ basis set, TIP4P-Ew~\\cite{Horn2004a} with point charges representing ions, for MBIS-S~\\cite{Verstraelen2016a}, for the CHARMM drude model of water (SWM4-NDP~\\cite{Lamoureux2006a}) with ions due to Yu {\\em et al.}~\\cite{Yu2010a}, as well as point core+Gaussian vsite (PC+GV), and point charge + Gaussian shell (PC+GS) derived here using ACT."
    label = "tab:ion_water2"
    wtable(file_path, models, mydata, caption, label)
    return file_path
    
def water_ions_induction():
    mymp  = "../AlexandriaFF/MP2-aug-cc-pvtz.xml"
    models = { "CHARMM": { "ff": "../ForceFields/CharmmDrude.xml", "mp": mymp },
               "PC+GS": { "ff": "../AlexandriaFF/PC+GS-elec.xml", "mp": mymp } }
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
    caption = "Water-ion induction energies at distances close to their energy minimum. Distance r (\\AA) between ions and water oxygen/hydrogen from Experiment (ref.~\\citenum{Heyrovska2006a}), and minimized water dimer (ref.~\\citenum{temelso2011benchmark}). Induction energies are reported in kJ/mol from the SAPT2+(CCD)-$\\delta$MP2 method with an aug-cc-pVTZ basis set, for the CHARMM drude model of water (SWM4-NDP~\\cite{Lamoureux2006a}) with ions due to Yu {\\em et al.}~\\cite{Yu2010a}, as well the point charge + Gaussian shell (PC+GS) derived here using ACT."
    label = "tab:ion_inducs"
    wtable(file_path, models, mydata, caption, label)
    return file_path
    
def ah_ions():
    models = { "PC": { "ff": "../ForceFields/TIP4PEW-JC.xml", "mp": mymp },
               "MBIS-S": { "ff": "../AlexandriaFF/P+S_updated.xml", "mp": None },
               "Walz":   { "ff": "../ForceFields/Walz2018a.xml", "mp": mymp },
               "PC+GV": { "ff": "../AlexandriaFF/PC+GV-elec.xml", "mp": mymp },
               "PC+GS": { "ff": "../AlexandriaFF/PC+GS-elec.xml", "mp": mymp } }
    newfn = "data-ah-ions.json"
    if not os.path.exists(newfn):
        sys.exit("Cannot find %s" % newfn)
    with open(newfn, "r") as inf:
        mydata = json.load(inf)

    add_calcs(mydata, models)

    file_path = "Ions-sapt2-JC-Walz2018a-ACT.tex"
    caption= "Ion-pair electrostatic energies at distances close to their energy minimum. Distance r (\\AA) between ions and electrostatic energies from the SAPT2+(CCD)$\\delta$MP2/aug-cc-pVTZ level of theory, for point charges (PC), for MBIS-S~\\cite{Verstraelen2016a}, for the Walz {\\em et al.} model with a Gaussian charge distribution~\\cite{Walz2018a}, and the ACT models PC+GV and PC+GS (see Methods). The RMSD and MSE were calculated with respect to the SAPT2+(CCD)$\\delta$MP2 with the aug-cc-pVTZ basis set electrostatic energy. Note that this level of theory is different from Tables~S1 and~S2 and the results cannot be compared directly. In addition, the results in Tables~S1 and~S2 are from fitting models to the ESP, whereas in this table training was done on SAPT data as indicated above."

    label = "tab:ion_ah"
    wtable(file_path, models, mydata, caption, label)
    return file_path
    

def ac_mt_gaff():
    models = { "RESP": { "ff": "../AlexandriaFF/coul-p.xml", "mp": mymp, "qtype": "qRESP" },
               "BCC": { "ff": "../AlexandriaFF/coul-p.xml", "mp": mymp, "qtype": "qBCC" },
               "MBIS-S": { "ff": "../AlexandriaFF/P+S_updated.xml", "mp": "../AlexandriaFF/sapt-esp5-mbiss.xml"  },
               "PC+GV": { "ff": "../AlexandriaFF/PC+GV-elec.xml", "mp": mymp },
               "PC+GS": { "ff": "../AlexandriaFF/PC+GS-elec.xml", "mp": mymp } }
    newfn = "data-sc-ions.json"
    if not os.path.exists(newfn):
        sys.exit("Cannot find %s" % newfn)
    with open(newfn, "r") as inf:
        mydata = json.load(inf)
    add_calcs(mydata, models)

    file_path = "AC-MA-IONS-GAFF.tex"
    caption = "Electrostatic energy (kJ/mol) between alkali ions, halides or water (oxygen) and amino acid side chain analogs, formate (oxygen), acetate (oxygen), methylammonium (nitrogen), ethylammonium (nitrogen) from SAPT2+(CCD)$\\delta$MP2/aug-cc-pVTZ, and charges determined using either RESP~\\cite{Bayly1993a} or BCC~\\cite{Jakalian2000a} as well as two models generated using the ACT."
    label = "tab:ac_ma_ions"
    wtable(file_path, models, mydata, caption, label)
    return file_path
    
if __name__ == "__main__":
    files = []
    #files.append(water_ions())
    #files.append(water_ions_induction())
    files.append(ac_mt_gaff())
    #files.append(ah_ions())
    print("Please check files:")
    for fn in files:
        print("  %s" % fn)
