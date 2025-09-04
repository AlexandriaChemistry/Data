#!/usr/bin/env python3

import os, glob, sys

def run_command(command):
    os.system(command)

train   = "Train"
test    = "Test"

acmparm = {
    "Mulliken":  { "ref": "Mulliken1955a", "ff": "coul-p.xml" },
    "Hirshfeld": { "ref": "Hirshfeld1977a",  "ff": "coul-p.xml" },
    "ESP":       { "ref": "Besler1990a", "ff": "coul-p.xml" },
    "CM5":       { "ref": "Marenich2012a",  "ff": "coul-p.xml" },
    "BCC":       { "ref": "Jakalian2000a", "ff": "coul-p.xml" },
    "RESP":      { "ref": "Bayly1993a", "ff": "coul-p.xml" },
    "ACM-esp-G":      { "ff": "esp-g.xml", "nparm": 48, "label": "GC", "target": "ESP" },
    "ACM-esp-GV":     { "ff": "esp-gv.xml", "nparm": 54, "label": "GC+PGV", "target": "ESP" },
    "ACM-elec-P":     { "ff": "coul-p.xml", "nparm": 32, "label": "PC", "target": "Elec" },
    "ACM-allelec-P":  { "ff": "all-p.xml", "nparm": 32, "label": "PC", "target": "Elec+Induc" },
    "ACM-elec-G":     { "ff": "coul-g.xml", "nparm": 48, "label": "GC", "target": "Elec" },
    "ACM-allelec-G":  { "ff": "all-g.xml", "nparm": 48, "label": "GC", "target": "Elec+Induc" },
    "ACM-elec-GV":    { "ff": "coul-gv.xml", "nparm": 62, "label": "GC+PGV", "target": "Elec" },
    "ACM-allelec-GV": { "ff": "all-gv.xml", "nparm": 62, "label": "GC+PGV", "target": "Elec+Induc" },
    "ACM-all-PG":     { "ff": "all-pg.xml", "nparm": 125, "label": "PC+GVS", "target": "Elec,Induc" }
}

def run_one(qtype:str) -> dict:
    if not qtype in acmparm:
        sys.exit("Unknown qtype %s" % qtype)
    molprops = "../AlexandriaFF/sapt-esp.xml"

    log_filename = f"{qtype}.log"
    base_command = f"alexandria train_ff -nooptimize -g {log_filename} -sel ../Selection/ac-total.dat -mp {molprops} -ff ../AlexandriaFF/{acmparm[qtype]['ff']}"

    print(f"Running command for {qtype}")
    if "ACM" in qtype:
        mycmd = base_command + " -charges ../AlexandriaFF/hf-aug-cc-pvtz.xml "
    else:
        mycmd = base_command + f" -qtype q{qtype} -charges ../AlexandriaFF/esp-paper-gaussian.xml "
    run_command(mycmd)

    print(f"Reading log file {log_filename}")
    mydict = {}
    for dataset in [ train, test ]:
        mydict[dataset] = {"COUL": {}, "ALLELEC": {}}
    with open(log_filename, "r") as inf:
        dset = None
        for line in inf:
            if "Results for" in line:
                words = line.strip().split()
                dset = words[3]
            elif "COULOMB (kJ" in line:
                words = line.strip().split()
                try:
                    mydict[dset]["COUL"]["N"] = int(words[2])
                    mydict[dset]["COUL"]["RMSD"] = round(float(words[5]),1)
                    mydict[dset]["COUL"]["MSE"] = round(float(words[6]),1)
                except ValueError:
                    sys.exit(f"Strange line {line.strip()}")
            elif "ALLELEC (kJ" in line:
                words = line.strip().split()
                try:
                    mydict[dset]["ALLELEC"]["N"] = int(words[2])
                    mydict[dset]["ALLELEC"]["RMSD"] = round(float(words[5]),1)
                    mydict[dset]["ALLELEC"]["MSE"] = round(float(words[6]),1)
                except ValueError:
                    sys.exit(f"Strange line {line.strip()}")

    return mydict

def get_train_test(logfn:str):
    ntrain = 0
    ntest  = 0
    with open(logfn, "r") as inf:
        for line in inf:
            if "COULOMB (kJ/mol)" in line:
                words = line.split()
                if "Train" in line:
                    ntrain = int(words[2])
                else:
                    ntest = int(words[2])
    return ntrain, ntest

mytable = {}
couls = ""
allelecs = ""
labels = ""

charge_models = [
    ("header", "Existing charge models" ),
    ("Mulliken", ""),
    ("Hirshfeld", ""),
    ("ESP", ""),
    ("CM5", ""),
    ("BCC", ""),
    ("RESP",""),
    ("header", "Non-polarizable ACT models based on monomer ESP" ),
    ("ACM", "-esp-G"), ("ACM", "-esp-GV"),
    ("header", "Non-polarizable ACT models based on SAPT dimer energies" ),
    ("ACM", "-elec-P"), ("ACM", "-allelec-P"),
    ("ACM", "-elec-G"), ("ACM", "-allelec-G"),
    ("ACM", "-elec-GV"), ("ACM", "-allelec-GV"),
    ("header", "Polarizable SAPT-based ACT model" ),
    ("ACM", "-all-PG")
]

for qt, suffix in charge_models:
    if qt == "header":
        continue
    qtsuf = qt+suffix
    mytable[qtsuf] = run_one(qtsuf)
    newcoul = f"COULOMB-{qtsuf}.xvg"
    run_command(f"mv COULOMB.xvg {newcoul}")
    couls += f" {newcoul}"
    newallelec = f"ALLELEC-{qtsuf}.xvg"
    run_command(f"mv ALLELEC.xvg {newallelec}")
    allelecs += f" {newallelec}"
    labels += f" {qtsuf}"
    for fn in [ "EXCHANGE.xvg", "EXCHIND.xvg", "DISPERSION.xvg", "INDUCTIONCORRECTION.xvg", "EPOT.xvg", "INDUCTION.xvg" ]:
        if os.path.exists(fn):
            os.unlink(fn)

run_command(f"viewxvg -f {couls} -label {labels} -ls None -mk o + x -res -noshow -pdf legacy_coul.pdf")
run_command(f"viewxvg -f {allelecs} -label {labels} -ls None -mk o + x -res -noshow -pdf legacy_allelec.pdf")

ntrain, ntest = get_train_test("ESP.log")

with open("legacy.tex", "w") as outf:
    outf.write("\\begin{table}[htb]\n")
    outf.write("\\centering\n")
    outf.write("\\caption{Root mean square deviation (RMSE) and mean signed error (MSE) of electrostatic energies (Elec, kJ/mol) and the sum of electrostatics and induction (Elec+Induc, kJ/mol) for popular charge models compared to SAPT2+(CCD)$\\delta$MP2 with the aug-cc-pVTZ basis set. The dataset consisted of 77 dimers (Table S6). \\#P indicates the number of parameters in the model. The training targets are indicated and RMSD and MSE values that correspond to the training set are indicated in {\\bf bold font}. A non-polarizable model with virtual sites with a Gaussian distributed charge (on anions and potassium ion only) is labeled as GC+PGV. The polarizable point charge + Gaussian virtual site and shell (PC+GVS) model was trained on electrostatic and induction energies in one step.}\n")
    
    
    outf.write("\\label{legacy}\n")
    outf.write("\\begin{tabular}{lcccccccc}\n")
    outf.write("\\hline\n")
    outf.write(" & Dataset & Training & \\#P & \\multicolumn{2}{c}{Elec}  & \\multicolumn{2}{c}{Elec+Induc}\\\\\n")
    outf.write("Model & & target & & RMSD & MSE & RMSD & MSE \\\\\n")

    for qt, suffix in charge_models:
        if qt == "header":
            outf.write("\\hline\n")
            outf.write("\\multicolumn{8}{c}{\\bf %s}\\\\\n" % suffix)
            continue
        qtsuf = qt + suffix
        label = qt
        if "label" in acmparm[qtsuf]:
            label = acmparm[qtsuf]["label"]
        star  = None
        if "target" in acmparm[qtsuf]:
            star  = acmparm[qtsuf]["target"]

        rmsd     = 'RMSD'
        mse      = 'MSE'
        na       = 'N/A'
        rmsd_str = {}
        mse_str  = {}
        for mydata in [ "COUL", "ALLELEC" ]:
            rmsd_str[mydata] = {}
            mse_str[mydata]  = {}
            for dataset in [ train, test ]:
                ttable = mytable[qtsuf][dataset][mydata]
                if rmsd in ttable and mse in ttable:
                    bold    = False
                    if star and "nparm" in acmparm[qtsuf] and dataset == train:
                        if "Elec,Induc" in star:
                            bold = True
                        elif "Elec+Induc" in star:
                            bold = mydata == "ALLELEC"
                        else:
                            bold = mydata == "COUL"
                        
                    if  bold:
                        rmsd_str[mydata][dataset] = f"\\textbf{{{ttable[rmsd]}}}"
                        mse_str[mydata][dataset]  = f"\\textbf{{{ttable[mse]}}}"
                    else:
                        rmsd_str[mydata][dataset] = f"{ttable[rmsd]}"
                        mse_str[mydata][dataset]  = f"{ttable[mse]}"
                else:
                    print("Something wrong with table for %s" % qtsuf)
                    sys.exit(ttable)
                np = ""
                if "nparm" in acmparm[qtsuf] and dataset == train:
                    np = acmparm[qtsuf]["nparm"]
        N = mytable[qtsuf][train][mydata]["N"]
        print(rmsd_str)
        print(mse_str)
        for dataset in [ train, test ]:
            target = ""
            if star:
                target = star
            cite = ""
            if "ref" in acmparm[qtsuf] and dataset == train:
                cite = f"~\\cite{{{acmparm[qtsuf]['ref']}}}"
            outf.write(f"{label}{cite} & {dataset} &{target} & {np} & {rmsd_str['COUL'][dataset]} & {mse_str['COUL'][dataset]} & {rmsd_str['ALLELEC'][dataset]} & {mse_str['ALLELEC'][dataset]} \\\\\n")

    outf.write("\\hline\n")
    outf.write("\\end{tabular}\n")
    outf.write("\\end{table}\n")
