#!/usr/bin/env python3

import json, math, os, sys

prefix    = "../../ACTdata/Training/ions-nobles-gases/TT2b/allhh/pol/vsite/"
viewflags = "-lfs 32 -alfs 32 -tickfs 28"
sapt      = "../../Psi4_ACT/sapt2+3(ccd)dmp2-aug-cc-pvtz/dimer-scans/"
debug     = False

def read_log_rmsd(filenm:str, ff:str)->dict:
    mylist = { "Test": {}, "Train": {} }
    toread = { "Elec": [ "COULOMB", "INDUCTION" ],
               "Induc": [ "INDUCTION", "DISPERSION", "EXCHANGE", "EPOT" ] }

    if os.path.exists(filenm):
        with open(filenm, "r") as inf:
            for line in inf:
                for ener in toread[ff]:
                    if line.startswith(f"{ener} (kJ/mol)"):
                        words = line.strip().split()
                        rmsd = float(words[5])
                        dset = words[8].split("-")[1]
                        mylist[dset][ener] = { "RMSD": rmsd, "N": int(words[2]) }

    return mylist

def runit(fn:str, logfn:str, ff:str, sel:int):
    actdata = "../../ACTdata/"
    sel = actdata + f"Training/ions-nobles-gases{sel}.dat"
    os.system(f"alexandria train_ff -ff {fn}/Train-ff_{ff}.xml -mp {actdata}MolProps/sapt2+3-final2 -sel {sel} -nooptimize -g {logfn}")
    
def run_stats(mk:str, logdir:str):
    os.makedirs(logdir, exist_ok=True)
    indpol = "Induction (polarization only)"
    indall = "Induction (complete)"
    rename = { "COULOMB": "Electrostatics", "EPOT": "Total" }
    mylist = {}
    allsel = [1, 2, 3]
    for sel in allsel:
        for repl in [ "A" ]:
            tdir = f"{prefix}TT2b-pg-{mk}-{sel}-{repl}"
            mylist[tdir] = {}
            for ff in [ "Elec", "Induc" ]:
                logfn = f"{logdir}/{ff}-{sel}-{repl}.log"
                runit(tdir, logfn, ff, sel)
                for xvg in [ "COULOMB", "INDUCTION", "DISPERSION", "EXCHANGE", "EPOT" ]:
                    fff = xvg + ".xvg"
                    if os.path.exists(fff):
                        os.system(f"mv {fff} {logdir}/{ff}-{xvg}-{sel}-{repl}.xvg")
                for xvg in [ "EXCHIND.xvg", "ALLELEC.xvg" ]:
                    os.unlink(xvg)
                if os.path.exists(logfn):
                    allset = read_log_rmsd(logfn, ff)
                    for dset in allset:
                        if not dset in mylist[tdir]:
                            mylist[tdir][dset] = {}
                        for ener in allset[dset]:
                            newener = ener
                            if ener in rename:
                                newener = rename[ener]
                            elif ener == "INDUCTION":
                                if ff == "Elec":
                                    newener = indpol
                                else:
                                    newener = indall
                            else:
                                newener = ener[0] + ener[1:].lower()
                            mylist[tdir][dset][newener] = allset[dset][ener]
    texfn = "stats.tex"
    print(mylist)
    with open(texfn, "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{Root mean square deviation from SAPT energies (kJ/mol) for three different selections (Table S1) separated in Train and Test sets.}\n")
        outf.write("\\label{tab:stats}\n")
        outf.write("\\begin{tabular}{lcccc}\n")
        outf.write("\\hline\n")
        outf.write("Term & & Set 1 & Set 2 & Set 3\\\\\n")
        outf.write("\\hline\n")
        order = [ "Electrostatics", indpol, indall,
                  "Dispersion", "Exchange", "Total" ]
        for ener in order:
            for dset in [ "Train", "Test" ]:
                if "Train" == dset:
                    outf.write("%s & %s" % ( ener, dset ))
                else:
                    outf.write(" & %s" % ( dset ))
                for fn in mylist:
                    if ener in mylist[fn][dset]:
                        outf.write(" & %.1f" % ( mylist[fn][dset][ener]["RMSD"]) )
                    else:
                        outf.write(" & ")
                outf.write("\\\\\n")
        outf.write("\\hline\n")
        outf.write("\\end{tabular}\n")
        outf.write("\\end{table}\n")
    print("Please check %s" % texfn)

def mydistance(dimer:str, ww:str)->float:
    resdir = sapt + dimer + "/" + ww
    if not os.path.exists(resdir):
        sys.exit("Cannot find dir '%s'" % resdir)
    js     = resdir + "/results.json"
    if not os.path.exists(js):
        print("No such file %s" % js)
        return 0
    with open(js, "r") as inf:
        mydata = json.load(inf)
    mols = "mols"
    if len(mydata[mols]) != 2:
        sys.exit("Incomprehensible number of molecules %d in %s" % ( len(mydata[mols]), js) )
    atoms  = "atoms"
    coords = "coords"
    # Squared mindist
    md2   = 1e8
    for n1 in mydata[mols][0][atoms]:
        for n2 in mydata[mols][1][atoms]:
            mmm2 = 0
            for m in range(3):
                mmm2 += (n1[coords][m]-n2[coords][m])**2
            md2 = min(md2, mmm2)
    return math.sqrt(md2)
    
def read_log_ener(filenm:str, mylist:dict)->dict:
    print("Will read %s" % filenm)
    moldata = {}
    with open(filenm, "r") as inf:
        found = False
        dataset = None
        dimer = None
        for line in inf:
            if line.startswith("Molecule "):
                words = line.strip().split()
                dataset = words[13]
                dimer = words[3][:-1]
                moldata[dimer] = []
                if debug:
                    print("Found dimer '%s' dataset %s" % ( dimer, dataset ))
            elif line.find("QM       ACT        QM       ACT") >= 0:
                found = True
            elif line.find("EPOT RMSD") >= 0:
                found = False
                dataset = None
            elif found:
                words = line.strip().split()
                if len(words) >= 16:
                    if dataset:
                        if not words[4] == "x":
                            myind = mydistance(dimer, words[16])
                            if myind > 0:
                                mylist[dataset].append( ( float(words[4]), float(words[5]) ) )
                                moldata[dimer].append( ( myind, float(words[4]), float(words[5]) ) )
                        else:
                            if debug:
                                print("Dimer '%s' has incorrect values" % dimer)
    return mylist, moldata

def do_mols(logdir:str, logfn:str, fntype:str, resonly:bool):
    allsel = [ 2 ]
    xvgdir = "xvgs-" + logfn + "-" + fntype
    os.makedirs(xvgdir, exist_ok=True)
    mylist = { "Train": [], "Test": [] }
    for sel in allsel:
        for repl in [ "A" ]:
            logfn = logdir + f"/{logfn}-{sel}-{repl}.log"
            _, moldata = read_log_ener(logfn, mylist)
            for md in moldata.keys():
                xvgfn = xvgdir + "/" + md.replace("#", "-") + ".xvg"
                with open(xvgfn, "w") as outf:
                    outf.write("@ xaxis label \"Distance\"\n")
                    outf.write("@ yaxis label \"Induction (kJ/mole)\"\n")
                    setind = 0
                    if not resonly:
                        outf.write("@ s%d legend \"SAPT\"\n" % setind)
                        setind += 1
                        outf.write("@ s%d legend \"ACT\"\n" % setind)
                        setind += 1
                    outf.write("@ s%d legend \"Residual\"\n" % setind)
                    for dd in sorted(moldata[md]):
                        if resonly:
                            outf.write("%10g  %10g\n" % ( dd[0], dd[2]-dd[1] ) )
                        else:
                            outf.write("%10g  %10g  %10g  %10g\n" % ( dd[0], dd[1], dd[2], dd[2]-dd[1] ) )

def do_fig1(logdir:str, fntype:str):
    allsel = [1, 2, 3]
    mylist = { "Train": [], "Test": [] }
    for sel in allsel:
        for repl in [ "A" ]:
            tdir  = f"{prefix}TT2b-pg-{fntype}-{sel}-{repl}"
            logfn = tdir + "/train-Elec.log"
            mylist, _ = read_log_ener(logfn, mylist)
            NN     = { "Train": 0, "Test": 0 }
            for ml in mylist.keys():
                xvgfn = f"induction{ml}.xvg"
                with open(xvgfn, "w") as outf:
                    outf.write("@ xaxis label \"Induction (kJ/mol)\"\n")
                    outf.write("@ yaxis label \"Alexandria\"\n")
                    for pair in mylist[ml]:
                        outf.write("%10g  %10g\n" % ( pair[0], pair[1] ) )
                        NN[ml] += 1
                print("Please check %s" % xvgfn)
                ltrain = f"Train (N = {NN['Train']})"
                ltest  = f"Test (N = {NN['Test']})"
            if sel == 2:
                elpdf = "electrostatics2.pdf"
                os.system(f"viewxvg -f {logdir}/Elec-COULOMB-2-A.xvg -ls None None -mk o x -res {viewflags} -pdf {elpdf}  -noshow -legend_y 0.3 -legend_x 0.65 ")
                print("Please check %s" % elpdf)
            indpdf = "induction2.pdf"
            os.system(f"viewxvg -f inductionTrain.xvg inductionTest.xvg -ls None None -mk o x  -label \"{ltrain}\" \"{ltest}\" -res -ymin -500 -ymax 500 -pdf {indpdf} -noshow -xmin -150 -legend_y 0.3 {viewflags}")

if __name__ == "__main__":
    mk     = "MORSEic-Kronecker3"
    logdir = "logs-" + mk
    run_stats(mk, logdir)
    do_fig1(logdir, mk)
    do_mols(logdir, "Elec", mk, True)
    do_mols(logdir, "Induc", mk, True)
    do_mols(logdir, "Induc", "MSic", True)
