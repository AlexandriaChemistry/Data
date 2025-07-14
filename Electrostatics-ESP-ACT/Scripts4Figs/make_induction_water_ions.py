#!/usr/bin/env python3

import os, sys

def run_one(ff:str, inducxvg:str, idlog:str):
    os.system("alexandria train_ff -nooptimize -ff %s -sel ../Selection/water+ions.dat -mp ../AlexandriaFF/sapt-esp.xml -g %s -charges ../AlexandriaFF/hf-aug-cc-pvtz" % ( ff, idlog ) )
    os.system("grep -v Alexandria INDUCTION.xvg | grep -v Train > %s" % inducxvg)
    
def rmsd(logfn:str)->float:
    with open(logfn, "r") as inf:
        for line in inf:
            if "INDUCTION (kJ/mol)" in line:
                words = line.strip().split()
                try:
                    return float(words[5])
                except ValueError:
                    continue
    sys.exit("Incomprehensible file %s" % logfn)
        
if __name__ == "__main__":
    ffs = [ { "ff": "../ForceFields/CharmmDrude.xml",
              "abbrev": "CharmmDrude", "label": "CHARMM Drude" },
              { "ff": "../AlexandriaFF/all-pg.xml",
              "abbrev": "PC+GVS", "label": "PC+GVS" }
           ]
    for ff in range(len(ffs)):
        logfn = ffs[ff]["abbrev"] + ".log"
        xvgfn = ffs[ff]["abbrev"] + ".xvg"
        run_one(ffs[ff]["ff"], xvgfn, logfn)
        ffs[ff]["xvg"] = xvgfn
        ffs[ff]["label"] += ( " RMSD %.1f kJ/mol" % ( rmsd(logfn) ) )

    pdf = "fig3.pdf"
    os.system("viewxvg -ls None -mk o x -alfs 32 -lfs 32 -tickfs 28 -f %s %s -label '%s' '%s' -pdf %s -noshow -res -legend_x 0.26 -legend_y 0.26" % ( ffs[0]["xvg"], ffs[1]["xvg"], ffs[0]["label"], ffs[1]["label"], pdf ) )
    print("Please check %s" % pdf)
