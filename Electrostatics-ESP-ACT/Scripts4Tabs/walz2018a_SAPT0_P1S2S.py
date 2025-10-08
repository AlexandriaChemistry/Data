#!/usr/bin/env python3

import json, math
import numpy as np
from potential_elec_functions import Point_core_1slater_2slater_shell, Point_core_1slater_shell, Point_core_2gaussian_shell, one_4pi_eps0, Point_core_gaussian_shell, Point_core_1slater_shell

sapt = {
    "Li-F":  { "rmin": 1.564, "eelec":  -867.15},
    "Li-Cl": { "rmin": 2.021, "eelec":  -661.43},
    "Li-Br": { "rmin": 2.17,  "eelec":  -612.35},
    "Na-F":  { "rmin": 1.926, "eelec":  -748.97},
    "Na-Cl": { "rmin": 2.361, "eelec":  -607.68},
    "Na-Br": { "rmin": 2.502, "eelec":  -572.36},
    "K-F":   { "rmin": 2.171, "eelec":  -704.64},
    "K-Cl":  { "rmin": 2.667, "eelec":  -566.07},
    "K-Br":  { "rmin": 2.821, "eelec":  -535.34}, }

walz = {
    "Li": { "zeta": 1.3772, "qcore": 3.391 },
    "Na": { "zeta": 1.2609, "qcore": 4.007 },
    "K":  { "zeta": 1.2412, "qcore": 5.464 },
    "F":  { "zeta": 0.9139, "qcore": 5.590 },
    "Cl": { "zeta": 0.6967, "qcore": 5.713 },
    "Br": { "zeta": 0.6465, "qcore": 5.778 } }

def gauss(distance:float, z1:float, z2:float)->float:
    zeta = z1*z1/math.sqrt(z1*z1+z2*z2)
    return one_4pi_eps0*math.erf(zeta*distance)/distance


def msd(ref_values, values):
    return (np.mean((np.array(ref_values) - np.array(values))**2))


def mse(ref_values, values):
    return np.mean(np.subtract(values,ref_values))

def compute_one(params:dict, a1:str, a2:str, rrr:float, func:int)->float:
    sf = str(func)
    if func in [ 0, 1]:
        if func == 0:
            myfunc = Point_core_gaussian_shell
        elif func == 1:
            myfunc = Point_core_1slater_shell
        [eesp] = myfunc([rrr],
                        params[a1]["q_c_"+sf], params[a1]["q_s_"+sf],
                        params[a2]["q_c_"+sf], params[a2]["q_s_"+sf],
                        params[a1]["z2_"+sf], params[a2]["z2_"+sf])
    elif func in [ 2, 3]:
        if func == 3:
            myfunc = Point_core_1slater_2slater_shell
        elif func == 2:
            myfunc = Point_core_2gaussian_shell
        [eesp]  = myfunc([rrr],
                         params[a1]["q_c_"+sf], params[a1]["q_s1_"+sf], params[a1]["q_s2_"+sf],
                         params[a2]["q_c_"+sf], params[a2]["q_s1_"+sf], params[a2]["q_s2_"+sf],
                         params[a1]["z1_"+sf], params[a2]["z1_"+sf],
                         params[a1]["z2_"+sf], params[a2]["z2_"+sf])
    else:
        sys.exit("Cannot handle func_index %d" % func)
    return eesp
        
def print_walz(outf, print_funcs:list, T:int):
    with open('../AnalyticalFitting/params_4_%s.json' % T, 'r') as json_f:
        params = json.load(json_f)

    epc_val  = { "msd": [], "mse": [] }
    eesp_val = {}
    for pf in print_funcs:
        eesp_val[pf] = { "msd": [], "mse": [] }

    for rm in sapt.keys():
        ions = rm.split("-")
        a1   = ions[0]
        a2   = ions[1]
        rrr  = sapt[rm]["rmin"]
        qc1  = walz[a1]["qcore"]
        qs1  = 1-qc1
        qc2  = walz[a2]["qcore"]
        qs2  = -1-qc2
        eelec = -gauss(rrr, walz[a1]["zeta"], walz[a2]["zeta"])

        # Point charge first
        epc     = -one_4pi_eps0/rrr
        diff    = epc - sapt[rm]["eelec"]
        epc_val["msd"].append(diff**2)
        epc_val["mse"].append(diff)
        outf.write("%s & %g & %.1f & %.1f " % ( rm.replace("-",""), rrr, sapt[rm]["eelec"], epc ) )
        # Then our own functions
        for pf in print_funcs:
            eesp  = compute_one(params, a1, a2, rrr, pf)
            diff  = eesp - sapt[rm]["eelec"]
            eesp_val[pf]["msd"].append(diff**2)
            eesp_val[pf]["mse"].append(diff)
            outf.write(" & %.1f" % eesp)

        outf.write("\\\\\n")

    outf.write("\\hline\n")
    outf.write("%s & & & %.1f" % ( "RMSD", np.sqrt(np.mean(epc_val["msd"])) ) )
    for pf in print_funcs:
        outf.write(" & %.1f" % np.sqrt(np.mean(eesp_val[pf]["msd"])) )
    outf.write("\\\\\n")
    outf.write("%s & & & %.1f" % ( "MSE", np.mean(epc_val["mse"]) ) )
    for pf in print_funcs:
        outf.write(" & %.1f" % np.mean(eesp_val[pf]["mse"]) )
    outf.write("\\\\\n")

if __name__ == "__main__":
    print_funcs = [ 0, 1, 2, 3 ]
    cols        = "c" * len(print_funcs)
    for T in [ 10, 100 ]:
        filename = ( "electable-%d.tex" % T)

        ctext       = { 0: "  a Gaussian charge", 1: "  a 1S Slater charge", 2: "  two Gaussian charges", 3: "  a 1S and a 2S Slater charge" }
        caption     = "Electrostatic energies at the experimental minimum energy distance~\\cite{NISTa} based on the SAPT0~\\cite{Sherill2013a} level of theory. Point charge (PC) energy follows from Coulomb's law. ESP indicates model consisting of a point charge, combined with "
    
        abcd   = "ABCD"
        myesps = ""
        for ppp in range(len(print_funcs)):
            pf = print_funcs[ppp]
            myesps += (" & ESP %c " % abcd[ppp] )
            caption += " " + str(abcd[ppp]) + ":"
            caption += ctext[pf]
            if ppp < len(print_funcs)-2:
                caption += ","
            elif ppp == len(print_funcs)-2:
                caption += " or "
        if T == 10:
            caption += " fitted to the Hartree-Fock electrostatic potential from 2.0 to 4.5 {\\AA}."
        elif T == 100:
            caption += " fitted to the Hartree-Fock electrostatic potential from 0.0 to 4.5 {\\AA}."
        with open(filename, "w") as outf:
            outf.write("""\\begin{table}[ht]
\\centering
\\caption{%s}
\\label{tab:sapt%d}
\\begin{tabular}{lccc%s}
\\hline
Ion pair & r$_{min}$ & \\multicolumn{%d}{c}{E$_{elec}$ (kJ/mol)}\\\\
& ({\\AA})&                 SAPT0 & PC %s\\\\
\\hline
""" % ( caption, T, cols, 2+len(print_funcs), myesps ) )
            print_walz(outf, print_funcs, T)
            outf.write("""\\hline
\\end{tabular}
\\end{table}
""")

        print("Output written to %s" % filename)
