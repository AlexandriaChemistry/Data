#!/usr/bin/env python3

import json, math, os, sys
from scipy.optimize import curve_fit
import numpy as np

prefix    = "../../ACTdata/Training/ions-nobles-gases/TT2b/allhh/pol/vsite/"
viewflags = "-lfs 32 -alfs 32 -tickfs 28"
sapt      = "../../Psi4_ACT/sapt2+3(ccd)dmp2-aug-cc-pvtz/dimer-scans/"
debug     = False

def morse(x, De, beta, bondlength):
    bx = beta*(x-bondlength)
    return De*( ( 1 - np.exp(-bx))**2 - 1)

def yukawa(x, A, B, C, D):
    return A*np.exp(-B*x)/x + C*np.exp(-D*x)

def fit_residual(dimer:str, distances, residual)->dict:
    cfit = { "Morse": { "func": morse,
                        "p0": [ 10, 1, 3 ],
                        "bounds": ( [ 0, 0.1, 0 ], [ 1000, 8, 10 ] ),
                        "params": None,
                        "value": None,
                        "msd": 0,
                        "rmsd": 0 },
             "Yukawa": { "func": yukawa,
                         "p0": [ 10, 1, -10, 1 ],
                         "bounds": ( [ -1000, 0.1, -100, 0.1 ], [ 1000, 5, 100, 5 ] ),
                         "params": None,
                         "value": None,
                         "msd": 0,
                         "rmsd": 0 }
            }
    print("There are %d data points for %s" % ( len(distances), dimer ) )
    pots = list(cfit.keys())
    for pot in pots:
        if debug:
            xx = 0.2
            print(f"morse_fit({xx}, {p0Morse}) = %g" % ( morse(xx, *p0Morse)))
        try:
            cfit[pot]["params"], _ = curve_fit(cfit[pot]["func"], distances, residual,
                                               p0=cfit[pot]["p0"], 
                                               bounds=cfit[pot]["bounds"],
                                               maxfev=10000)
            cfit[pot]["values"] = []
            index = 0
            for d in distances:
                val = cfit[pot]["func"](d, *(cfit[pot]["params"]))
                cfit[pot]["values"].append(val)
                cfit[pot]["msd"] += (residual[index]-val)**2
                index += 1
      
        except Exception as e:
            print(f"Error fitting {pot} data for {dimer}: {e}")

        cfit[pot]["rmsd"] = math.sqrt(cfit[pot]["msd"]/len(distances))

    return cfit

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

def print_morse_parameters(pot:str, mfit:dict, texfn:str):
    with open(texfn, "w") as outf:
        outf.write("\\begin{longtable}{lcccc}\n")
        outf.write("\\caption{Parameters for fit of induction residual to a %s potential}\\\\\n" % pot)
        outf.write("\\hline\n")
        outf.write("Dimer & De & $\\beta$ & r$_{min}$ & RMSD (kJ/mole)\\\\\n")
        param = "param"
        for md in sorted(mfit.keys()):
            outf.write("%s " % md.replace("#", "-") )
            for p in mfit[md][param]:
                outf.write(" & %.2f " % p)
            outf.write(" %.1f \\\\\n" % ( mfit[md]["rmsd"] ) )
        outf.write("\\hline\n")
        outf.write("\\end{longtable}\n")
        
def do_mols(logdir:str, logfn:str, fntype:str, resonly:bool, dofit:bool):
    allsel = [ 2 ]
    xvgdir = "xvgs-" + logfn + "-" + fntype
    os.makedirs(xvgdir, exist_ok=True)
    mylist = { "Train": [], "Test": [] }
    if dofit:
        mfit = {}
    for sel in allsel:
        for repl in [ "A" ]:
            logfn = logdir + f"/{logfn}-{sel}-{repl}.log"
            _, moldata = read_log_ener(logfn, mylist)
            for md in moldata.keys():
                xvgfn = xvgdir + "/" + md.replace("#", "-") + ".xvg"
                # Prepare data for fitting
                if dofit:
                    distances = []
                    residual  = []
                    for dd in sorted(moldata[md]):
                        distances.append(dd[0])
                        # The residual is SAPT - Polarization, that is the missing part of the energy function
                        residual.append( dd[1]-dd[2] )
                    cfit = fit_residual(md, distances, residual)
                    for pot in cfit:
                        if not pot in mfit:
                            mfit[pot] = {}
                        if "params" in cfit[pot] and not cfit[pot]["params"] is None:
                            mfit[pot][md] = { "param": cfit[pot]["params"].copy(), "rmsd": cfit[pot]["rmsd"] }
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
                    if dofit:
                        for pot in mfit:
                            if not "values" in cfit[pot]:
                                continue
                            setind += 1
                            outf.write("@ s%d legend \"%s fit to residual, RMSD = %.1f\"\n" %
                                       ( setind, pot, mfit[pot][md]["rmsd"] ) )

                    if dofit:
                        for index in range(len(distances)):
                            outf.write("%10g  %10g" % ( distances[index], residual[index] ) )
                            for pot in mfit:
                                if "values" in cfit[pot]:
                                    outf.write("  %10g" % cfit[pot]["values"][index] )
                            outf.write("\n")
                    else:
                        for dd in sorted(moldata[md]):
                            if resonly:
                                outf.write("%10g  %10g  %10g\n" % ( dd[0], dd[1], dd[2] ) )
                            else:
                                outf.write("%10g  %10g  %10g  %10g\n" % ( dd[0], dd[1], dd[2], dd[2]-dd[1] ) )
    if dofit:
        for pot in mfit:
            print_morse_parameters(pot, mfit[pot], f"{pot}_param.tex")

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
#    run_stats(mk, logdir)
#    do_fig1(logdir, mk)
    do_mols(logdir, "Elec", mk, True, True)
    do_mols(logdir, "Induc", mk, True, False)
    do_mols(logdir, "Induc", "MSic", True, False)
