#!/usr/bin/env python3

import json, math
from Electrostatic_4models import Point_core_1slater_2slater_shell

one_4pi_eps0 = 1389

sapt = {
    "Li-F":  { "rmin": 1.564, "eelec": -850.53137313 },
    "Li-Cl": { "rmin": 2.021, "eelec": -656.84841560 },
    "Li-Br": { "rmin": 2.17, "eelec": -609.54342356 },
    "Na-F":  { "rmin":1.926, "eelec": -750.61771962 },
    "Na-Cl": { "rmin": 2.361, "eelec": -608.0 },
    "Na-Br": { "rmin": 2.502, "eelec": -572.68092126 },
    "K-F":   { "rmin":2.171, "eelec": -719.80224985 },
    "K-Cl":  { "rmin": 2.667, "eelec": -570.18777426 },
    "K-Br":  { "rmin": 2.821, "eelec": -538.37985053 }, }
    
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
    
def print_walz(outf):
    with open('../AnalyticalFitting/params_4.json', 'r') as json_f:
        params = json.load(json_f)

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
        # q_c_na, q_s1_na, q_s2_na, q_c_cl, q_s1_cl, q_s2_cl, z_c_na, z_c_cl, z_s_na, z_s_cl):
        [eesp]  = Point_core_1slater_2slater_shell([rrr], 
                                                   params[a1]["q_c_"+a1+"_3"], params[a1]["q_s1_"+a1+"_3"], params[a1]["q_s2_"+a1+"_3"],
                                                   params[a2]["q_c_"+a2+"_3"], params[a2]["q_s1_"+a2+"_3"], params[a2]["q_s2_"+a2+"_3"],
                                                   params[a1]["z1_"+a1+"_3"], params[a1]["z2_"+a1+"_3"],
                                                   params[a2]["z1_"+a2+"_3"], params[a2]["z2_"+a2+"_3"])
        epc   = -one_4pi_eps0/rrr
        outf.write("%s & %g & %.0f & %.0f & %.0f & %.0f\\\\\n" % ( rm, rrr, sapt[rm]["eelec"], epc, eelec, eesp) )

if __name__ == "__main__":
    filename = "electable.tex"
    with open(filename, "w") as outf:
        outf.write("""\\begin{table}[ht]
        \\centering
        \\begin{tabular}{lccccc}
        \\hline
        Ion pair & r$_{min}$ & \\multicolumn{4}{c}{E$_{elec}$ (kJ/mol)}\\\\
        & ({\\AA})&                 SAPT & PC & Walz & ESP \\\\
        \\hline
        """)
        print_walz(outf)
        outf.write("""\\hline
    \\end{tabular}
\\caption{Electrostatic energies from SAPT at the experimental minimum energy distance~\\cite{NISTa}  and corresponding energy of two point charges. The 
SAPT2+(CCD)$\\delta$MP2~\\cite{Sherill2013a} level of theory was applied and the total Electrostatics energy is reported. Point charge (PC) energy follows 
from Coulomb's law. Walz indicates a previous model with a Gaussian charge distribution~\\cite{Walz2018a}. ESP indicates model consisting of a point charge, a 1S Slater and a 2S Slater distribution, fitted to the electrostatic potential from 0.1 to 4.5 {\\AA}.}
    \\label{tab:sapt}
    \\end{table}
    """)

    print("Output written to %s" % filename)
