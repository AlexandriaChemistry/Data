#!/usr/bin/env python3

import os, sys

mpdef   = " -mp ../AlexandriaFF/MP2-aug-cc-pvtz.xml"
qdef    = " -charges ../AlexandriaFF/MP2-aug-cc-pvtz.xml"
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
    { "label": "PC+GV3*", "ff": "PC+GV-esp3", "flags": " -qalg ESP -mp ../AlexandriaFF/MP2-aug-cc-pvtz_Updated.xml", train: monomer },
    { "label": "PC+SV3*", "ff": "PC+SV-esp3", "flags": " -qalg ESP -mp ../AlexandriaFF/MP2-aug-cc-pvtz_Updated.xml", train: monomer },
    { "label": "PC+GV4*", "ff": "PC+GV-esp4", "flags": " -qalg ESP -mp ../AlexandriaFF/MP2-aug-cc-pvtz_Updated.xml", train: monomer },
    { "label": "PC+SV4*", "ff": "PC+SV-esp4", "flags": " -qalg ESP -mp ../AlexandriaFF/MP2-aug-cc-pvtz_Updated.xml", train: monomer },
    { "label": "PC", "ff": "PC-elec", "flags":  mpdef, train: dimer },
    { "label": "GC", "ff": "GC-elec", "flags":  mpdef, train: dimer },
    { "label": "SC", "ff": "SC-elec", "flags":  mpdef, train: dimer },
    { "label": "PC+GV", "ff": "PC+GV-elec", "flags":  mpdef, train: dimer },
    { "label": "PC+SV", "ff": "PC+SV-elec", "flags":  mpdef, train: dimer },
    { "label": "PC+GS", "ff": "PC+GS-elec", "flags":  mpdef, train: dimer }
]

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
    return rmsd, mylog
    
def get_mol_rmsd(mylog:str)->dict:
    rmsd = {}
    with open(mylog, "r") as inf:
        lines = inf.readlines()
        for ii in range(len(lines)):
            line = lines[ii].strip()
            
            if line.startswith("Molecule "):
                words = line.split()
                if len(words) < 14 or ii+2 >= len(lines):
                    continue
                mol = words[3][:-1]
                try:
                    if lines[ii+2].startswith("ESP rms"):
                        www = lines[ii+2].strip().split()
                        rmsd[mol] = float(www[2])
                except ValueError:
                    print("Cannot read ESP RMSD for %s in %s" % ( mol, mylog ))
    return rmsd

def write_table(tab:str, myrmsd:dict):
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

def write_table_mol(tab:str, rmsd_mol:dict):
    mols = []
    for model in rmsd_mol:
        if len(rmsd_mol[model].keys()) > len(mols):
            mols = list(rmsd_mol[model].keys())

    with open(tab, "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{RMSD (kJ/mol e) with respect to the electrostatic potential computed at the MP2/aug-cc-pvtz level of theory (see Methods) for the side chain analogs in Table 1 and water. For description of the different models and training, see Methods.}\n")
        outf.write("\\label{tab:esprms}\n")
        outf.write("\\begin{tabular}{lc")
        for c in range(1+len(mols)):
            outf.write("c")
        outf.write("}\n")
        outf.write("\\hline\n")
        outf.write("Model & Train ")
        for c in range(len(mols)):
            outf.write(" & \\rotatebox{90}{%s}" % mols[c])
        outf.write("&\\rotatebox{90}{Average}\\\\\n")
        outf.write("\\hline\n")
        for model in rmsd_mol:
            ttt = None
            for k in range(len(allffs)):
                if allffs[k]["label"] == model:
                    ttt = allffs[k][train]
                    break
            outf.write("%s & %s " % ( model, ttt ))
            rmsdtot = 0
            nrmsd   = 0
            for mol in mols:
                if mol in rmsd_mol[model]:
                    outf.write("& %.1f" % ( rmsd_mol[model][mol] ))
                    rmsdtot += rmsd_mol[model][mol]
                    nrmsd   += 1
                else:
                    outf.write(" & ")
            outf.write(" & %.1f \\\\\n" % ( rmsdtot / nrmsd ) )
        outf.write("\\hline\n")
        outf.write("\\end{tabular}\n")
        outf.write("\\end{table}\n")
    print("Please check %s" % tab)

if __name__ == "__main__":
    myrmsd = { monomer: [], dimer: [] }
    rmsd_mol = {}
    for qqm in allffs:
        print("Will run %s" % qqm['label'])
        ff = f"../AlexandriaFF/{qqm['ff']}.xml"
        rmsd, mylog   = run_one(ff, qqm['label'], qqm["flags"])
        rmsd_mol[qqm["label"]] = get_mol_rmsd(mylog)
        if rmsd:
            myrmsd[qqm[train]].append( ( qqm['label'],  rmsd ) )

    for tt in myrmsd:
        for kv in myrmsd[tt]:
            print("%s  %g" % ( kv[0], kv[1] ))
    write_table("esprmsd.tex", myrmsd)
    write_table_mol("esprmsd_mol.tex", rmsd_mol)
