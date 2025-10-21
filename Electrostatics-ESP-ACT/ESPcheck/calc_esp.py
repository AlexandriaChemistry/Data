#!/usr/bin/env python3

import os, sys

mpdef = " -mp ../AlexandriaFF/MP2-aug-cc-pvtz.xml"
qdef  = " -charges ../AlexandriaFF/MP2-aug-cc-pvtz.xml"

def run_one(ff:str, myid:str, flags:str)->float:
    myff  = os.path.basename(ff)
    mylog = f"{myid}.log"
    cmd = ("alexandria train_ff -ff %s -nooptimize -sel espmono2.dat -g '%s' -v 4 %s" % ( ff, mylog, flags ) )
    print("Command = '%s'" % cmd)
    sys.stdout.flush()
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
        { "label": "Mulliken", "ff": "coul-p", "flags": " -qqm qMulliken "+mpdef+qdef, train: monomer },
        { "label": "Hirshfeld", "ff": "coul-p", "flags": " -qqm qHirshfeld "+mpdef+qdef, train: monomer },
        { "label": "ESP", "ff": "coul-p", "flags": " -qqm qESP "+mpdef+qdef, train: monomer },
        { "label": "CM5", "ff": "coul-p", "flags":  " -qqm qCM5 "+mpdef+qdef, train: monomer },
        { "label": "BCC", "ff": "coul-p", "flags":  " -qqm  qBCC "+mpdef+qdef, train: monomer },
        { "label": "RESP", "ff": "coul-p", "flags":  " -qqm qRESP "+mpdef+qdef, train: monomer },
        { "label": "MBIS", "ff": "coul-p", "flags":  " -qqm qMBIS  -charges ../AlexandriaFF/MP2-MBIS.xml"+mpdef, train: monomer },
        { "label": "MBIS-S", "ff": "P+S_updated", "flags": " -qalg None -mp ../AlexandriaFF/MP2-aug-cc-pvtz_Updated.xml",  train: monomer },
        { "label": "PC+GV*", "ff": "PC+GV-esp", "flags": " -qalg ESP -mp ../AlexandriaFF/MP2-aug-cc-pvtz_Updated.xml", train: monomer },
        { "label": "PC+SV*", "ff": "PC+SV-esp", "flags": " -qalg ESP -mp ../AlexandriaFF/MP2-aug-cc-pvtz_Updated.xml", train: monomer },
        { "label": "PC", "ff": "PC-elec", "flags":  mpdef, train: dimer },
        { "label": "GC", "ff": "GC-elec", "flags":  mpdef, train: dimer },
        { "label": "SC", "ff": "SC-elec", "flags":  mpdef, train: dimer },
        { "label": "PC+GV", "ff": "PC+GV-elec", "flags":  mpdef, train: dimer },
        { "label": "PC+SV", "ff": "PC+SV-elec", "flags":  mpdef, train: dimer },
        { "label": "PC+GS", "ff": "PC+GS-elec", "flags":  mpdef, train: dimer }
    ]
    myrmsd = { monomer: [], dimer: [] }
    for qqm in allffs:
        print("Will run %s" % qqm['label'])
        ff = f"../AlexandriaFF/{qqm['ff']}.xml"
        rmsd = run_one(ff, qqm['label'], qqm["flags"])

        if rmsd:
            myrmsd[qqm[train]].append( ( qqm['label'],  rmsd ) )

    for tt in myrmsd:
        for kv in myrmsd[tt]:
            print("%s  %g" % ( kv[0], kv[1] ))
    tab = "esprmsd.tex"
    with open(tab, "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{RMSD (kJ/mol e) of different charge models with respect to the electrostatic potential computed at the MP2/aug-cc-pvtz level of theory (see Methods) for the side chain analogs in Table 1 and water. For description of the different models and training, see Methods.}\n")
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
