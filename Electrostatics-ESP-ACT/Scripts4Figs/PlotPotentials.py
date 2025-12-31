#!/usr/bin/env python3

import math, os

BOHR  = 0.527 # Angstrom
KJMOL = 1389

def gen_pot(xi:float, r0:float, r1:float, n:int, qtot:int, func:str):
    mypot = []
    qc    = 1 
    for i in range(n):
        r = r0 + i*(r1-r0)/(n-1)
        if func == "Slater":
            y = qc/r + (qtot-qc)*(1/r-(1+r*xi)*math.exp(-2*r*xi)/(r))
        elif func == "Gaussian":
            y = qc/r + (qtot-qc)*math.erf(r*xi)/(r)
        elif func == "PC":
            y = qtot/r
        mypot.append( ( r, y * KJMOL ) )
    return mypot
    
if __name__ == "__main__":
    qtot = -1
    xi   = 0.8
    for func in [ "Slater", "Gaussian", "PC" ]:
        pot = gen_pot(xi, 0.5, 5, 81, qtot, func)
        xvg = ( "pot2-%s.xvg" % func )
        with open(xvg, "w") as outf:
            outf.write("@ xaxis label 'Distance (Angstrom)'\n")
            outf.write("@ yaxis label 'Potential (kJ/mol e)'\n")
            for p in pot:
                outf.write("%10f  %10f\n" % ( p[0], p[1] ) )
    pdf = "esp.pdf"
    os.system("plotxvg -f pot2-PC.xvg pot2-Gaussian.xvg pot2-Slater.xvg  -label PC PC+GC PC+SC  -alfs 28 -lfs 28 -tickfs 24 -legend_x 0.5 -legend_y 0.5 -pdf %s -noshow -mk" % pdf)
    print("Please check %s" % pdf)
