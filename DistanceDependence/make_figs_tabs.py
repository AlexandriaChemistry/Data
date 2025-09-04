#!/usr/bin/env python3

import json, math, os, sys, argparse
from scipy.optimize import curve_fit
import numpy as np

prefix    = "../../ACTdata/Training/ions-nobles-gases/TT2b/allhh/pol/vsite/"
viewflags = "-lfs 32 -alfs 32 -tickfs 28"
sapt      = "../../Psi4_ACT/sapt2+3(ccd)dmp2-aug-cc-pvtz/dimer-scans/"
debug     = False
Induction = "Induction"
Electrostatics = "Electrostatics"

def ms(x, A, b):
    return -A*np.exp(-b*x)

def morse(x, De, beta, bondlength):
    bx = beta*(x-bondlength)
    return De*( ( 1 - np.exp(-bx))**2 - 1)

def yukawa(x, A, B, C, D):
    return A*np.exp(-B*x)/x + C*np.exp(-D*x)

def fit_residual(dimer:str, distances, residual, init_p0)->dict:
    cfit = { 
        "MS": { "func": ms,
                "p0": [ 10, 1 ],
                "bounds": ( [ 0, 0.1 ], [ 5000, 12 ] ),
                "params": None,
                "value": None,
                "msd": 0,
                "rmsd": 0 },
        "Morse": { "func": morse,
                   "p0": [ 10, 1, 3 ],
                   "bounds": ( [ 0, 0.1, 0 ], [ 5000, 12, 20 ] ),
                   "params": None,
                   "value": None,
                   "msd": 0,
                   "rmsd": 0 },
        "Yukawa": { "func": yukawa,
                    "p0": [ 10, 1, -10, 1 ],
                    "bounds": ( [ -500000, 0.1, 0, 0.1 ], [ 0, 5, 500000, 5 ] ),
                    "params": None,
                    "value": None,
                    "msd": 0,
                    "rmsd": 0 }
    }
    if debug:
        print("There are %d data points for %s" % ( len(distances), dimer ) )
    pots = list(cfit.keys())
    for pot in pots:
        if debug:
            xx = 0.2
            print(f"morse_fit({xx}, {p0Morse}) = %g" % ( morse(xx, *p0Morse)))
        try:
            # Initial params may be passed from the calling routine
            init_param = cfit[pot]["p0"]
            if init_p0 and pot in init_p0:
                init_param = init_p0[pot]["param"]
            # Check bounds
            for p in range(len(init_param)):
                init_param[p] = max(init_param[p], cfit[pot]["bounds"][0][p])
                init_param[p] = min(init_param[p], cfit[pot]["bounds"][1][p])
            cfit[pot]["params"], _ = curve_fit(cfit[pot]["func"], distances, residual,
                                               p0=init_param,
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
               "Elec+Induc": [  "COULOMB", "INDUCTION" ],
               "Inter": [ "DISPERSION", "EXCHANGE", "EPOT" ] }

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

def read_log_rmsd_aver(filenm:str)->dict:
    mylist = { "Train": {}, "Test": {} }

    with open(filenm, "r") as inf:
        dataset = None
        dimer = None
        for line in inf:
            if line.startswith("Molecule "):
                words = line.strip().split()
                dataset = words[13]
                dimer = words[3][:-1]
                if debug:
                    print("Found dimer '%s' dataset %s" % ( dimer, dataset ))
            elif line.find(" RMSD") >= 0 and dimer and dataset:
                words = line.strip().split()
                if len(words) >= 14:
                    try:
                        ener = words[0]
                        if not ener in mylist[dataset]:
                            mylist[dataset][ener] = { "RMSD": 0, "N": 0 }
                        mylist[dataset][ener]["RMSD"] += float(words[2])
                        mylist[dataset][ener]["N"]    += 1
                    except ValueError:
                        print("Cannot interpret line '%s'"  % (line.strip()))
            elif dimer and len(line)==0:
                dimer = None
                dataset = None

    allset = { "Train": {}, "Test": {} }
    for dset in mylist:
        for ener in mylist[dset]:
            if ener == "N":
                continue
            allset[dset][ener] = { "RMSD": mylist[dset][ener]["RMSD"]/mylist[dset][ener]["N"],
                                   "N": mylist[dset][ener]["N"] }
    return allset

def read_log_ener(filenm:str, induc:bool)->dict:
    mylist = { "Train": [], "Test": [] }
    print("Will read %s" % filenm)
    moldata = {}
    # Which energies to fetch for induction
    iqm = 4
    iact = 5
    if not induc:
        # Get the electrostatics
        iqm = 2
        iact = 3
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
                        if not words[iqm] == "x":
                            myind = mydistance(dimer, words[16])
                            if myind > 0:
                                mylist[dataset].append( ( float(words[iqm]), float(words[iact]) ) )
                                moldata[dimer].append( ( myind, float(words[iqm]), float(words[iact]) ) )
                        else:
                            if debug:
                                print("Dimer '%s' has incorrect values" % dimer)
    return mylist, moldata


def run_stats(mk:str, logdir:str, allsel:list, repls: list)->int:
    def run_train_ff(fn:str, logfn:str, ff:str, sel:int):
        actdata = "../../ACTdata/"
        sel = actdata + f"Training/ions-nobles-gases{sel}.dat"
        xml = f"{fn}/Train-ff_{ff}.xml"
        if os.path.exists(xml):
            return os.system(f"alexandria train_ff -ff {xml} -mp {actdata}MolProps/sapt2+3-final2 -sel {sel} -nooptimize -g {logfn}")
        return 1

    os.makedirs(logdir, exist_ok=True)
    indpol = "Induction\\footnote{Polarization only}"
    indall = "Induction\\footnote{Including correction}"
    rename = { "COULOMB": Electrostatics, "EPOT": "Total" }
    mylist = {}
    aver_rmsd = True
    for sel in allsel:
        for repl in repls:
            tdir = f"{prefix}TT2b-pg-{mk}-{sel}-{repl}"
            mylist[tdir] = {}
            for ff in [ "Elec", "Elec+Induc", "Inter" ]:
                logfn = f"{logdir}/{ff}-{sel}-{repl}.log"
                if run_train_ff(tdir, logfn, ff, sel) != 0:
                    continue
                for xvg in [ "COULOMB", "INDUCTION", "DISPERSION", "EXCHANGE", "EPOT" ]:
                    fff = xvg + ".xvg"
                    if os.path.exists(fff):
                        os.system(f"mv {fff} {logdir}/{ff}-{xvg}-{sel}-{repl}.xvg")
                for xvg in [ "EXCHIND.xvg", "ALLELEC.xvg" ]:
                    if os.path.exists(xvg):
                        os.unlink(xvg)
                if os.path.exists(logfn):
                    if aver_rmsd:
                        allset = read_log_rmsd_aver(logfn)
                    else:
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
        outf.write("\\caption{Root mean square deviation from SAPT energies (kJ/mol) for five different selections (Table S1) separated in Train and Test sets.}\n")
        outf.write("\\label{tab:stats}\n")
        outf.write("\\begin{tabular}{lc")
        for r in allsel:
            outf.write("c")
        outf.write("}\n")
        outf.write("\\hline\n")
        outf.write("Term &")
        for r in range(len(allsel)):
            outf.write("& Set %d" % (r+1))
        outf.write("\\\\\n")
        outf.write("\\hline\n")
        order = [ Electrostatics, indpol, indall,
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
    
def print_parameters(mfit:dict):
    invaa = "(1/{\\AA})"
    kjm   = "(kJ/mole)"
    rmin  = "r$_{min}$ ({\\AA})"
    texprops = {  "MS": { "columns": "lccc", "heading": 
                             f"Dimer & A {kjm} & $\\beta$ {invaa} & RMSD {kjm}" },
                  "Morse": { "columns": "lcccc", "heading": 
                             f"Dimer & De {kjm} & $\\beta$ {invaa} & {rmin} & RMSD {kjm}" },
                  "Yukawa": { "columns": "lccccc", "heading":
                              f"Dimer & A {kjm}& B {invaa} & C {kjm} & D {invaa} & RMSD {kjm}" }
                }
    for pot in texprops.keys():
        texfn = pot + "_param.tex"
        with open(texfn, "w") as outf:
            outf.write("\\begin{longtable}{%s}\n" % texprops[pot]["columns"])
            outf.write("\\caption{Parameters for fit of a %s potential to the induction residual (SAPT - Polarization).}\\\\\n" % pot)
            outf.write("\\hline\n")
            outf.write("%s\\\\\n" % texprops[pot]["heading"])
            outf.write("\\hline\n")
            param = "param"
            msd    = 0
            nmsd   = 0
            nparam = 0
            for md in sorted(mfit.keys()):
                if pot in mfit[md]:
                    if nparam == 0:
                        nparam = len(mfit[md][pot][param])
                    outf.write("%s " % md.replace("#", "-") )
                    for p in mfit[md][pot][param]:
                        outf.write(" & %.2f " % p)
                    outf.write(" & %.2f \\\\\n" % ( mfit[md][pot]["rmsd"] ) )
                    msd += mfit[md][pot]["rmsd"]**2
                    nmsd += 1
            outf.write("\\hline\n")
            outf.write("Average")
            for p in range(nparam):
                outf.write(" & ")
            outf.write(" & %.2f\\\\\n\\hline\n" % ( math.sqrt(msd/nmsd) ) )
            outf.write("\\end{longtable}\n")

def make_pdf(xvgdir:str, plotfn:str, xvgfns:list):
    count = 0
    with open(plotfn, "w") as outf:
        for md in sorted(xvgfns.keys()):
            if len(xvgfns[md]) < 2:
                continue
            dimer = md.replace("#", "-")
            pdf   = f"{xvgdir}/{dimer}.pdf"
            # Plot both in a panel of two
            os.system(f"viewxvg -f {xvgfns[md][Electrostatics]}  {xvgfns[md][Induction]} -pdf {pdf} -xframe 9 -yframe 12 -legend_x 0.3 -legend_y 0.7 -panels -noshow")
            if os.path.exists(pdf):
                outf.write("\\begin{figure}[ht]\n")
                outf.write("\\centering\n")
                outf.write("\\includegraphics[width=14cm]{%s}\n" % pdf)
                outf.write("\\caption{A) residual electrostatic energy (SAPT - ACT) for dimer {\\bf %s}, B) difference between SAPT induction energy and analytical curve(s) fitted to it.}\n\\end{figure}\n\n" % ( dimer ) )
                count += 1
            if count == 10:
                outf.write("\\cleardoublepage\n")
                count = 0

def write_xvg(xvgfn:str, eterm:str, resonly:bool, args,
              ener_data:dict, mfit:dict, cfit:dict,
              distances, residual)->int:
    setind = 0
    with open(xvgfn, "w") as outf:
        outf.write("@ xaxis label \"Distance\"\n")
        outf.write("@ yaxis label \"%s (kJ/mole)\"\n" % eterm)
        if not resonly:
            outf.write("@ s%d legend \"SAPT\"\n" % setind)
            setind += 1
            outf.write("@ s%d legend \"ACT\"\n" % setind)
            setind += 1
        outf.write("@ s%d legend \"Residual\"\n" % setind)
        if args.fit:
            if None != mfit:
                for pot in mfit:
                    if not "values" in cfit[pot]:
                        continue
                    setind += 1
                    outf.write("@ s%d legend \"%s fit, RMSD = %.2f\"\n" %
                               ( setind, pot, mfit[pot]["rmsd"] ) )

            index = 0
            for dd in sorted(ener_data):
                outf.write("%10g  %10g" % ( dd[0], dd[1] - dd[2] ) )
                if None != mfit:
                    for pot in mfit:
                        if "values" in cfit[pot]:
                            outf.write("  %10g" % cfit[pot]["values"][index] )
                outf.write("\n")
                index += 1
        else:
            for dd in sorted(ener_data):
                if resonly:
                    outf.write("%10g  %10g  %10g\n" % ( dd[0], dd[1], dd[2] ) )
                else:
                    outf.write("%10g  %10g  %10g  %10g\n" % ( dd[0], dd[1], dd[2], dd[2]-dd[1] ) )
    return setind

def do_mols(logdir:str, logfn:str, fntype:str, resonly:bool, args, allsel:list,
            repls: list):
    xvgdir = "xvgs-" + logfn + "-" + fntype
    os.makedirs(xvgdir, exist_ok=True)
    mylist = { "Train": [], "Test": [] }
    if args.fit:
        mfit = {}
    if args.dist:
        dist_range = {}
    # Fetch initial fitting parameters if relevant
    initial_fit = None
    if args.load and args.fit:
        if len(args.load) > 0 and not os.path.exists(args.load):
            sys.exit("File '%s' does not exist" % args.load)
        with open(args.load, "r") as inf:
            initial_fit = json.load(inf)
    for sel in allsel:
        for repl in repls:
            mylogfn = logdir + f"/{logfn}-{sel}-{repl}.log"
            _, induc_data = read_log_ener(mylogfn, True)
            if args.plot:
                _, elec_data = read_log_ener(mylogfn, False)
            xvgfns = {}
            for md in induc_data.keys():
                # Prepare data for fitting
                if args.fit:
                    distances = []
                    residual  = []
                    mfit[md]  = {}
                    for dd in sorted(induc_data[md]):
                        distances.append(dd[0])
                        # The residual is SAPT - Polarization, that is the missing part of the energy function
                        residual.append( dd[1]-dd[2] )
                    # Check distance range for table S1
                    # Has to be based on elec_data
                    if args.dist and not md in dist_range:
                        mydd = sorted(elec_data[md])
                        dist_range[md] = ( mydd[0][0], mydd[-1][0] )
                    init_p0 = None
                    if initial_fit and md in initial_fit:
                        init_p0 = initial_fit[md]
                    cfit = fit_residual(md, distances, residual, init_p0)
                    for pot in cfit:
                        if "params" in cfit[pot] and not cfit[pot]["params"] is None:
                            mfit[md][pot] = { "param": cfit[pot]["params"].tolist(), "rmsd": cfit[pot]["rmsd"] }
                if args.plot:
                    xvgfns[md] = {}
                    dimer = md.replace("#", "-")
                    xvgfn1 = xvgdir + f"/{Induction}-{dimer}.xvg"
                    write_xvg(xvgfn1, Induction, resonly, args, induc_data[md],
                              mfit[md], cfit, distances, residual)
                    xvgfn2 = xvgdir + f"/{Electrostatics}-{dimer}.xvg"
                    write_xvg(xvgfn2, Electrostatics, True, args, elec_data[md],
                              None, None, None, None)
                    if os.path.exists(xvgfn1) and os.path.exists(xvgfn2):
                        xvgfns[md] = { Induction: xvgfn1, Electrostatics: xvgfn2 }
            if args.plot:
                make_pdf(xvgdir, args.plot, xvgfns)
    if args.fit:
        print_parameters(mfit)
        if args.save:
            with open(args.save, "w") as outf:
                json.dump(mfit, outf, sort_keys=True, indent=4)
                print("Stored parameters to %s" % args.save)

    if args.dist:
        with open(args.dist, "w") as outf:
            outf.write("#Dimer,MinDist,MaxDist\n")
            for md in sorted(dist_range.keys()):
                outf.write("%s,%g,%g\n" % ( md, dist_range[md][0], dist_range[md][1] ) )
            print("Stored distance range to %s" % args.dist)

def do_fig1(logdir:str, allsel:list, repls:list):
    for sel in allsel:
        for repl in repls:
            pdf = f"{logdir}/fig1-{sel}-{repl}.pdf"            
            os.system(f"viewxvg -f {logdir}/Elec-COULOMB-{sel}-{repl}.xvg {logdir}/Elec-INDUCTION-{sel}-{repl}.xvg  -pdf {pdf} -xframe 9 -yframe 12 -legend_x 0.6 -legend_y 0.3 -panels -noshow -res -ls None None -mk o + -xframe 10 -yframe 14  -lfs 26 -tickfs 20 -alfs 26 ")
            print("Please check %s" % pdf)

def parser():
    parser  = argparse.ArgumentParser(description="""
Generate figures and tables by doing analysis on ACT training outputs.
You have to run the -stats option (or -all) at least once before using the other flags.
""")
    parser.add_argument("-all", "--all", help="Do everything", action="store_true")
    parser.add_argument("-stats", "--stats", help="Run statistics and make table", action="store_true")
    parser.add_argument("-pol", "--pol", help="Analyse induction just based on polarization", action="store_true")
    parser.add_argument("-fit", "--fit", help="Fit residual (SAPT-ACT polarization) to analytical functions (activates the -pol flag as well)", action="store_true")
    defload = "params.json"
    parser.add_argument("-load", "--load", help="Load initial fitting parameters from, default "+defload, type=str, default=defload)
    parser.add_argument("-save", "--save", help="Save the fitting parameters to...", type=str, default=None)
    parser.add_argument("-induc", "--induc", help="Analyse induction just based on polarization", action="store_true")
    parser.add_argument("-dist", "--dist", help="Save distance range to this file", type=str,default=None)
    parser.add_argument("-plot", "--plot", help="Make pdfs of residual plots and produce latex include file with name ...", type=str, default=None)

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parser()
    
    mk     = "MORSEic-Kronecker-EXTRA"
    logdir = "logs-" + mk
    allsel = [ 1, 2, 3, 4, 5 ]
    repls  = [ "C" ]
    if args.stats or args.all:
        run_stats(mk, logdir, allsel, repls)
    if os.path.isdir(logdir):
        do_fig1(logdir, allsel, repls)
        if args.pol or args.all or args.fit:
            do_mols(logdir, "Elec", mk, True, args, allsel, repls)
        if args.induc or args.all:
            do_mols(logdir, "Elec+Induc", mk, True, args, allsel, repls)
            do_mols(logdir, "Elec+Induc", "MSic", True, args, allsel, repls)
