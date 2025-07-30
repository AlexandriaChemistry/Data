#!/usr/bin/env python3

import os, sys

def read_log(filenm:str)->dict:
    mylist = { "Test": [], "Train": [] }
    with open(filenm, "r") as inf:
        found = False
        dataset = None
        dimer = None
        for line in inf:
            if line.startswith("Molecule "):
                words = line.strip().split()
                dataset = words[13]
                dimer = words[3][:-1]
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
                            mylist[dataset].append( ( float(words[4]), float(words[5]) ) )
                        else:
                            print("Dimer '%s' has incorrect values" % dimer)
    return mylist

if __name__ == "__main__":
    viewflags = "-lfs 32 -alfs 32 -tickfs 28"
    prefix = "../../ACTdata/Training/ions-nobles-gases/TT2b/allhh/pol/vsite/"
    for fn in [ "800elec3", "800elec4", "800elec5" ]:
        logfn = prefix + fn + "/train-Elec.log"
        mylist = read_log(logfn)
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
        if fn == "800elec4":
            elpdf = "electrostatics4.pdf"
            os.system(f"viewxvg -f ../../ACTdata/Training/ions-nobles-gases/TT2b/allhh/pol/vsite/800elec4/COULOMB.xvg -ls None None -mk o x -res {viewflags} -pdf {elpdf}  -noshow -legend_y 0.3 -legend_x 0.65 ")
            print("Please check %s" % elpdf)
        os.system(f"viewxvg -f inductionTrain.xvg inductionTest.xvg -ls None None -mk o x  -label \"{ltrain}\" \"{ltest}\" -res -ymin -500 -ymax 500 -pdf {fn}.pdf -noshow -xmin -150 -legend_y 0.3 {viewflags}")
