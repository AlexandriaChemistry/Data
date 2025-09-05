#!/usr/bin/env python3

import os

def grep_rms(filter:str)->dict:
    filenm = ( "rmsdump-%s.txt" % filter )
    os.system("cd compounds; zgrep 'Charges from ESP fit,' */*-%s.log.gz > ../%s" % ( filter, filenm ) )

    mydict = {}
    with open(filenm, "r", encoding='utf-8') as inf:
        for line in inf:
            words = line.split("/")
            mol   = words[0]
            try:
                rms   = float(words[1].split()[6])
                if not mol in mydict or rms < mydict[mol]:
                    mydict[mol] = rms
            except ValueError:
                print("Could not read line '%s'" % line.strip())
    return mydict
    
if __name__ == "__main__":
    #filter = "B3LYP-aug-cc-pVTZ"
    filter = "HF-6-311G**"
    mydict = grep_rms(filter)
    maxrms = 0
    nmol   = 0
    for mol in mydict:
        rms = mydict[mol]
#        print("%s  %g" % ( mol, rms ))
        maxrms  = max(maxrms, rms)
        nmol   += 1
    mymax = 0.005
    print("Largest RMS found %g in %d compounds but using %g" % (maxrms, len(mydict.keys()),  mymax))
    histo = [0 for i in range(101)]
    outlier = open("outlier.txt", "w")
    rsum = 0
    for mol in mydict:
        index = int(100*mydict[mol]/mymax)
        rsum += mydict[mol]
        if index < 100:
            histo[index] += 1
        else:
            outlier.write("Outlier %s RMSD %g\n" % ( mol, mydict[mol] ))
            histo[100] += 1
    outlier.close()
    hartree = 2625.5
    rsum /= (len(mydict.keys()) - histo[100])
    print("Please check %d outliers in outlier.txt" % (histo[100]))
    print("Average RMSD %g kJ/mol e" % (hartree*rsum))
    histofn = f"histo-{filter}.xvg"
    with open(histofn, "w") as outf:
        outf.write("@ xaxis label \"RMSD (kJ/mol e)\"\n")
        outf.write("@ yaxis label \"(arbitrary units)\"\n")
        for i in range(len(histo)-1):
            outf.write("%10g  %10g\n" % ( hartree*i*mymax/100, histo[i]*1.0/nmol))

    pdf = "rmsd_histo.pdf"
    os.system("viewxvg -f %s -lfs 40 -alfs 40 -tickfs 32 -noshow -pdf %s" % ( histofn, pdf ) )
    print("Please check %s" % pdf )
