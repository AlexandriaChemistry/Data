#!/usr/bin/env python3

import os

def run_one(ff:str, myid:str, qtype:str, charges:str)->float:
    myff  = os.path.basename(ff)
    mylog = f"{myid}.log"
    cmd = ("alexandria train_ff -ff %s -mp ../AlexandriaFF/hf-aug-cc-pvtz.xml -nooptimize -sel ../Selection/monomer.dat -g %s -v 4" % ( ff, mylog ) )
    if qtype:
        if not charges:
            charges = "../AlexandriaFF/hf-aug-cc-pvtz.xml"
        cmd += (" -qtype q%s -charges %s" % ( qtype, charges ) )
    os.system(cmd)
    rmsd = None
    with open(mylog, "r") as inf:
        for line in inf:
            if line.find("ESP (kJ") >= 0:
                words = line.split()
                try:
                    rmsd = float(words[6])
                except ValueError:
                    print("Cannot read ESP RMSD in %s" % mylog)
    return rmsd
    
if __name__ == "__main__":
    allffs = [
        { "label": "Mulliken", "ff": "coul-p", "qtype": "Mulliken" },
        { "label": "Hirshfeld", "ff": "coul-p", "qtype": "Hirshfeld" },
        { "label": "ESP", "ff": "coul-p", "qtype": "ESP" },
        { "label": "CM5", "ff": "coul-p", "qtype": "CM5" },
        { "label": "BCC", "ff": "coul-p", "qtype": "BCC" },
        { "label": "RESP", "ff": "coul-p", "qtype": "RESP" },
        { "label": "MBIS", "ff": "coul-p", "qtype": "RESP", "charges": "../AlexandriaFF/mbis_ccsd.xml" },
        { "label": "MBIS-S", "ff": "P+S-hacked", "qtype": "RESP", "charges": "../AlexandriaFF/mbisS_ccsd.xml" },
        { "label": "PC", "ff": "PC-elec", "qtype": None },
        { "label": "GC", "ff": "GC-elec", "qtype": None },
        { "label": "SC", "ff": "SC-elec", "qtype": None },
        { "label": "PC+GV", "ff": "PC+GV-elec", "qtype": None },
        { "label": "PC+SV", "ff": "PC+SV-elec", "qtype": None },
        { "label": "PC+GS", "ff": "PC+GS-elec", "qtype": None }
    ]
    myrmsd = []
    for qtype in allffs:
        print("Will run %s" % qtype['label'])
        ff = f"../AlexandriaFF/{qtype['ff']}.xml"
        charges = None
        if "charges" in qtype:
            charges = qtype["charges"]
        rmsd = run_one(ff, qtype['label'], qtype['qtype'], charges)

        if rmsd:
            myrmsd.append( ( qtype['label'],  rmsd ) )
    for kv in myrmsd:
        print("%s  %g" % ( kv[0], kv[1] ))
    tab = "esprmsd.tex"
    with open(tab, "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{RMSD (kJ/mol e) of different charge models with respect to the electrostatic potential computed at the HF/aug-cc-pvtz level of theory (see Methods). For description of the different models and training, see Methods.}\n")
        outf.write("\\label{tab:esprms}\n")
        outf.write("\\begin{tabular}{lcc}\n")
        outf.write("\\hline\n")
        outf.write("Model & Training target & RMSD \\\\\n")
        outf.write("\\hline\n")
        for kv in myrmsd:
            train = "elec"
            outf.write("%s & %s & %.1f\\\\\n" % ( kv[0], train, kv[1] ))
        outf.write("\\hline\n")
        outf.write("\\end{tabular}\n")
        outf.write("\\end{table}\n")
    print("Please check %s" % tab)
