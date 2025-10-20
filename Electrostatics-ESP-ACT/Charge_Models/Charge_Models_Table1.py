#!/usr/bin/env python3

import os, glob, sys

debug   = False
train   = "Train"
test    = "Test"
header  = "header"

# Number of independent charges disregarding symmetry
nq = 45
acmparm = {
    "Mulliken":      { "ref": "Mulliken1955a", "nparm": nq, "ff": "coul-p.xml", header: "Existing charge models" },
    "Hirshfeld":     { "ref": "Hirshfeld1977a", "nparm": nq, "ff": "coul-p.xml" },
    "ESP":           { "ref": "Besler1990a", "nparm": nq, "ff": "coul-p.xml" },
    "CM5":           { "ref": "Marenich2012a",  "nparm": nq, "ff": "coul-p.xml" },
    "BCC":           { "ref": "Jakalian2000a", "nparm": nq, "ff": "coul-p.xml" },
    "RESP":          { "ref": "Bayly1993a", "nparm": nq, "ff": "coul-p.xml" },
    "MBIS":          { "ref": "Verstraelen2016a", "nparm": nq, "ff": "coul-p.xml" },
    "MBIS-S":        { "ref": "Verstraelen2016a", "nparm": nq*3, "ff": "P+S_updated.xml" },
    "PC+GV-esp":     { "ff": "PC+GV-esp.xml", "nparm": nq*3, "label": "PC+GV", "target": "ESP", header: "Non-polarizable ACT monomer-based models" },
    "PC+SV-esp":     { "ff": "PC+SV-esp.xml", "nparm": nq*3, "label": "PC+SV", "target": "ESP" },
    "PC-elec":       { "ff": "PC-elec.xml", "nparm": 66, "label": "PC", "target": "Elec", header: "Non-polarizable ACT dimer-based models" },
    "PC-allelec":    { "ff": "PC-allelec.xml", "nparm": 66, "label": "PC", "target": "Elec+Induc" },
    "GC-elec":       { "ff": "GC-elec.xml", "nparm": 85, "label": "GC", "target": "Elec" },
    "GC-allelec":    { "ff": "GC-allelec.xml", "nparm": 85, "label": "GC", "target": "Elec+Induc" },
    "SC-elec":       { "ff": "SC-elec.xml", "nparm": 85, "label": "SC", "target": "Elec" },
    "SC-allelec":    { "ff": "SC-allelec.xml", "nparm": 85, "label": "SC", "target": "Elec+Induc" },
    "PC+GV-elec":    { "ff": "PC+GV-elec.xml", "nparm": 106, "label": "PC+GV", "target": "Elec" },
    "PC+GV-allelec": { "ff": "PC+GV-allelec.xml", "nparm": 106, "label": "PC+GV", "target": "Elec+Induc" },
    "PC+SV-elec":    { "ff": "PC+SV-elec.xml", "nparm": 106, "label": "PC+SV", "target": "Elec" },
    "PC+SV-allelec": { "ff": "PC+SV-allelec.xml", "nparm": 106, "label": "PC+SV", "target": "Elec+Induc" },
    "PC+GS-elec":    { "ff": "PC+GS-elec.xml", "nparm": 156, "label": "PC+GS", "target": "Elec,Induc", header: "Polarizable ACT dimer-based models" },
    "PC+GS-allelec": { "ff": "PC+GS-allelec.xml", "nparm": 156, "label": "PC+GS", "target": "Elec+Induc" }
}

def myround(word:str)->float:
    fff = float(word)
    if abs(fff) >= 10:
        return round(fff, 0)
    else:
        return round(fff, 1)

def run_one(qtype:str, qm:str) -> dict:
    if not qtype in acmparm:
        sys.exit("Unknown qtype %s" % qtype)
    molprops = "../AlexandriaFF/sapt-esp5.xml"
    if qtype == "MBIS-S":
        molprops = "../AlexandriaFF/sapt-esp5-mbiss.xml"

    log_filename = f"{qtype}_{qm}.log"
    base_command = f"alexandria train_ff -nooptimize -g {log_filename} -sel ../Selection/ac-total.dat -mp {molprops} -ff ../AlexandriaFF/{acmparm[qtype]['ff']} "

    qfn = f"../AlexandriaFF/{qm}-aug-cc-pvtz.xml"
    print(f"Running command for {qtype}")
    if "elec" in qtype:
        mycmd = base_command + f" -charges {qfn} -qalg SQE"
    elif "esp" in qtype:
        mycmd = base_command + f" -charges {qfn} -qalg ESP"
    elif qtype == "MBIS":
        mycmd = base_command + f" -qqm qMBIS -charges ../AlexandriaFF/MP2-MBIS.xml -qalg Read "
    elif qtype == "MBIS-S":
        mycmd = base_command + " -qqm None " #+ f" -charges ../AlexandriaFF/mbisS_ccsd.xml -qqm qRESP "
    else:
        mycmd = base_command + f" -qqm q{qtype} -charges {qfn} -qalg Read "
    os.system(mycmd)

    if not os.path.exists(log_filename):
        return {}

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
                    mydict[dset]["COUL"]["RMSD"] = myround(words[5])
                    mydict[dset]["COUL"]["MSE"] = myround(words[6])
                except ValueError:
                    sys.exit(f"Strange line {line.strip()}")
            elif "ALLELEC (kJ" in line:
                words = line.strip().split()
                try:
                    mydict[dset]["ALLELEC"]["N"] = int(words[2])
                    mydict[dset]["ALLELEC"]["RMSD"] = myround(words[5])
                    mydict[dset]["ALLELEC"]["MSE"] = myround(words[6])
                except ValueError:
                    sys.exit(f"Strange line {line.strip()}")

    return mydict

def get_train_test(logfn:str):
    ntrain = 0
    ntest  = 0
    if os.path.exists(logfn):
        with open(logfn, "r") as inf:
            for line in inf:
                if "COULOMB (kJ/mol)" in line:
                    words = line.split()
                    if "Train" in line:
                        ntrain = int(words[2])
                    else:
                        ntest = int(words[2])
    return ntrain, ntest

def do_all(qm:str):
    mytable = {}
    couls = ""
    allelecs = ""
    labels = ""

    for qtsuf in acmparm:
        mytable[qtsuf] = run_one(qtsuf, qm)
        # Check whether we got any data
        if not mytable[qtsuf]:
            continue
        newcoul = f"COULOMB-{qtsuf}.xvg"
        os.system(f"mv COULOMB.xvg {newcoul}")
        couls += f" {newcoul}"
        newallelec = f"ALLELEC-{qtsuf}.xvg"
        os.system(f"mv ALLELEC.xvg {newallelec}")
        allelecs += f" {newallelec}"
        labels += f" {qtsuf}"
        for fn in [ "EXCHANGE.xvg", "EXCHIND.xvg", "DISPERSION.xvg", "INDUCTIONCORRECTION.xvg", "EPOT.xvg", "INDUCTION.xvg" ]:
            if os.path.exists(fn):
                os.unlink(fn)

    cmd1 = f"viewxvg -f {couls} -label {labels} -ls None -mk o + x -res -noshow -save legacy_coul_{qm}.pdf -panel"
    if debug:
        print(f"cmd1 = {cmd1}")
    os.system(cmd1)
    cmd2 = f"viewxvg -f {allelecs} -label {labels} -ls None -mk o + x -res -noshow -save legacy_allelec_{qm}.pdf -panel"
    if debug:
        print(f"cmd2 = {cmd2}")
    os.system(cmd2)

    ntrain, ntest = get_train_test("ESP_{qm}.log")

    with open(f"legacy_{qm}.tex", "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{Root mean square deviation (RMSD) and mean signed error (MSE), both in kJ/mol of electrostatic energies (Elec) and the sum of electrostatics and induction (Elec+Induc) for popular charge models compared to ACT models based on ESP or SAPT (Table~\\ref{tab:models}). The dataset consisted of 94 dimers (Table S6) and \\#P indicates the number of parameters in the model. Values corresponding to the training targets are indicated in {\\bf bold font}.}\n")
    
        outf.write("\\label{tab:legacy}\n")
        outf.write("\\begin{tabular}{lccccccccccc}\n")
        outf.write("\\hline\n")
        outf.write(" & Training & \\#P & \\multicolumn{2}{c}{Elec}  & \\multicolumn{2}{c}{Elec+Induc}& \\multicolumn{2}{c}{Elec}  & \\multicolumn{2}{c}{Elec+Induc}\\\\\n")
        outf.write("Model & target & & RMSD & MSE & RMSD & MSE & RMSD & MSE & RMSD & MSE \\\\\n")
        outf.write("& & & \\multicolumn{4}{c}{Train}& \\multicolumn{4}{c}{Test}\\\\\n")

        for qtsuf in acmparm:
            if header in acmparm[qtsuf]:
                outf.write("\\hline\n")
                outf.write("&&\\multicolumn{9}{c}{\\bf %s}\\\\\n" % acmparm[qtsuf][header])
            label = qtsuf
            if "label" in acmparm[qtsuf]:
                label = acmparm[qtsuf]["label"]
            star  = None
            if "target" in acmparm[qtsuf]:
                star  = acmparm[qtsuf]["target"]

            rmsd     = 'RMSD'
            mse      = 'MSE'
            na       = 'N/A'
            np       = ""
            if "nparm" in acmparm[qtsuf]:
                np = acmparm[qtsuf]["nparm"]
            rmsd_str = {}
            mse_str  = {}
            for mydata in [ "COUL", "ALLELEC" ]:
                rmsd_str[mydata] = {}
                mse_str[mydata]  = {}
                for dataset in [ train, test ]:
                    if not qtsuf in mytable:
                        print(f"No {qtsuf} in data")
                        continue
                    if not dataset in mytable[qtsuf]:
                        print(f"No {dataset} data in {qtsuf}")
                        continue
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
                            rmsd_str[mydata][dataset] = f"\\textbf{{{ttable[rmsd]:g}}}"
                            mse_str[mydata][dataset]  = f"\\textbf{{{ttable[mse]:g}}}"
                        else:
                            rmsd_str[mydata][dataset] = f"{ttable[rmsd]:g}"
                            mse_str[mydata][dataset]  = f"{ttable[mse]:g}"
                    else:
                        print("Something wrong with table for %s" % qtsuf)
                        sys.exit(ttable)
            if train in mytable[qtsuf]:
                N = mytable[qtsuf][train][mydata]["N"]
                if debug:
                    print(rmsd_str)
                    print(mse_str)
                target = ""
                if star:
                    target = star
                cite = ""
                if "ref" in acmparm[qtsuf] and dataset == train:
                    cite = f"~\\cite{{{acmparm[qtsuf]['ref']}}}"
                thisnp = np
                outf.write(f"{label}{cite} & {target} & {thisnp} & {rmsd_str['COUL'][train]} & {mse_str['COUL'][train]} & {rmsd_str['ALLELEC'][train]} & {mse_str['ALLELEC'][train]} & {rmsd_str['COUL'][test]} & {mse_str['COUL'][test]} & {rmsd_str['ALLELEC'][test]} & {mse_str['ALLELEC'][test]} \\\\\n")

        outf.write("\\hline\n")
        outf.write("\\end{tabular}\n")
        outf.write("\\end{table}\n")

    print("ntrain %d ntest %d" % ( ntrain, ntest ))

if __name__ == "__main__":
    #do_all("HF")
    do_all("MP2")
