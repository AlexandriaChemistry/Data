#!/usr/bin/env python3

import os

def run_one(ff:str, qtype:str)->float:
    myff = os.path.basename(ff)
    mylog = myff.replace(".xml", ".log")
    cmd = ("alexandria train_ff -ff ../AlexandriaFF/%s -mp ../../../ACTdata/MolProps/hf-aug-cc-pvtz.xml -nooptimize -sel ../Selection/monomer.dat -g %s -v 4" % ( myff, mylog ) )
    if qtype:
        cmd += (" -qtype q%s -charges ../AlexandriaFF/hf-aug-cc-pvtz.xml " % qtype)
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
        "Mulliken",
        "Hirshfeld",
        "ESP",
        "CM5",
        "BCC",
        "RESP",
        "esp-g.xml",
        "esp-gv.xml",
        "coul-p.xml", 
        "all-p.xml", 
        "coul-g.xml",
        "all-g.xml",
        "coul-gv.xml",
        "all-gv.xml",
    	"all-pg.xml",
    ]
    label = {
        "esp-g.xml": { "model": "GC", "train": "ESP" },
        "esp-gv.xml": { "model": "GC+PGV", "train": "ESP" },
        "coul-p.xml": { "model": "PC", "train": "SAPT (Elec)" },
        "all-p.xml": { "model": "PC", "train": "SAPT (Elec+Induc)" },
        "coul-g.xml": { "model": "GC", "train": "SAPT (Elec)" },
        "all-g.xml": { "model": "GC", "train": "SAPT (Elec+Induc)" },
        "coul-gv.xml": { "model": "GC+PGV", "train": "SAPT (Elec)" },
        "all-gv.xml": { "model": "GC+PGV", "train": "SAPT (Elec+Induc)" },
    	"all-pg.xml": { "model": "PC+GVS", "train": "SAPT (Elec,Induc)" }
    }
    myrmsd = []
    for ff in allffs:
        myff  = "coul-p.xml"
        qtype = ff
        myid  = qtype
        if ff.endswith(".xml"):
            myff  = ff
            qtype = None
            myid  = myff
        print("Will run %s" % myff)
        rmsd = run_one(myff, qtype)
        print("==================")
        if rmsd:
            myrmsd.append( ( myid,  rmsd ) )
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
            model = kv[0]
            train = ""
            if kv[0] in label:
                model = label[kv[0]]["model"]
                train = label[kv[0]]["train"]
            outf.write("%s & %s & %.1f\\\\\n" % ( model, train, kv[1] ))
        outf.write("\\hline\n")
        outf.write("\\end{tabular}\n")
        outf.write("\\end{table}\n")
    print("Please check %s" % tab)
