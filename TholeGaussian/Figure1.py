#!/usr/bin/env python3

import math

def thole(r:float, aiw:float)->float:
    return 1 - (1 + r / (2 * aiw)) * math.exp(- r / aiw)
    
def gauss(r:float, beta:float)->float:
    return math.erf(beta * r)
    
if __name__ == "__main__":
    mytab = [ { "alpha": 0.0005, "beta": 12 },
#              { "alpha": 0.001, "beta": 10 },
              { "alpha": 0.002, "beta": 8 } ]
    for mt in mytab:
        alpha = mt["alpha"]
        tw    = 2.6
        aiw   = ((alpha**2)**(1.0/6.0))/tw
        beta  = 2/(3*aiw*math.sqrt(math.pi))
#        beta  = mt["beta"]
        outfn = ( "tg%g.xvg" % beta )
        with open(outfn, "w") as outf:
            outf.write("@ xaxis label \"Distance (nm)\"\n")
            outf.write("@ s0 legend \"Thole, alpha = %g nm$^3$\"\n" % alpha)
            outf.write("@ s1 legend \"Gaussian, zeta = %g nm$^{-1}$\"\n" % beta)
            for ir in range(36):
                r = 0.01*ir
                outf.write("%10g  %10g  %10g\n" %
                           ( r, thole(r, aiw), gauss(r, beta) ) )
