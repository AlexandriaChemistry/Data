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
    monomer = "Monomer"
    dimer   = "Dimer"
    train   = "train"
    allffs = [
        { "label": "Mulliken", "ff": "coul-p", "qtype": "Mulliken", train: monomer },
        { "label": "Hirshfeld", "ff": "coul-p", "qtype": "Hirshfeld", train: monomer },
        { "label": "ESP", "ff": "coul-p", "qtype": "ESP", train: monomer },
        { "label": "CM5", "ff": "coul-p", "qtype": "CM5", train: monomer },
        { "label": "BCC", "ff": "coul-p", "qtype": "BCC", train: monomer },
        { "label": "RESP", "ff": "coul-p", "qtype": "RESP", train: monomer },
        { "label": "MBIS", "ff": "coul-p", "qtype": "RESP", "charges": "../AlexandriaFF/mbis_ccsd.xml", train: monomer },
        { "label": "MBIS-S", "ff": "P+S_updated", "qtype": None, train: monomer  },
        { "label": "PC", "ff": "PC-elec", "qtype": None, train: dimer },
        { "label": "GC", "ff": "GC-elec", "qtype": None, train: dimer },
        { "label": "SC", "ff": "SC-elec", "qtype": None, train: dimer },
        { "label": "PC+GV", "ff": "PC+GV-elec", "qtype": None, train: dimer },
        { "label": "PC+SV", "ff": "PC+SV-elec", "qtype": None, train: dimer },
        { "label": "PC+GS", "ff": "PC+GS-elec", "qtype": None, train: dimer }
    ]
    myrmsd = { monomer: [], dimer: [] }
    for qtype in allffs:
        print("Will run %s" % qtype['label'])
        ff = f"../AlexandriaFF/{qtype['ff']}.xml"
        charges = None
        if "charges" in qtype:
            charges = qtype["charges"]
        rmsd = run_one(ff, qtype['label'], qtype['qtype'], charges)

        if rmsd:
            myrmsd[qtype[train]].append( ( qtype['label'],  rmsd ) )

    for tt in myrmsd:
        for kv in myrmsd[tt]:
            print("%s  %g" % ( kv[0], kv[1] ))
    tab = "esprmsd.tex"
    with open(tab, "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{RMSD (kJ/mol e) of different charge models with respect to the electrostatic potential computed at the HF/aug-cc-pvtz level of theory (see Methods) for all side chain analogs in Table 1 and water. For description of the different models and training, see Methods.}\n")
        outf.write("\\label{tab:esprms}\n")
        outf.write("\\begin{tabular}{lclc}\n")
        outf.write("\\hline\n")
        outf.write("Model & RMSD & Model & RMSD  \\\\\n")
        outf.write("\\multicolumn{2}{c}{Monomer-based} & \\multicolumn{2}{c}{Dimer-based} \\\\\n")
        outf.write("\\hline\n")
        for c in range(len(myrmsd[monomer])):
            outf.write("%s & %.1f " % ( myrmsd[monomer][c][0], myrmsd[monomer][c][1] ))
            if c < len(myrmsd[dimer]):
                outf.write("& %s & %.1f\\\\\n" % ( myrmsd[dimer][c][0], myrmsd[dimer][c][1] ))
            else:
                outf.write(" & & \\\\\n")
        outf.write("\\hline\n")
        outf.write("\\end{tabular}\n")
        outf.write("\\end{table}\n")
    print("Please check %s" % tab)
