#!/usr/bin/env python3

import json, os

from walz2018a import one_4pi_eps0

with open('output_4_100.json', 'r') as json_f:
    output_data = json.load(json_f)

start = 50
xref  = output_data["data"]["Br"][0]
pc    = [0]*start
for i in range(start,len(xref)):
    pc.append(one_4pi_eps0/xref[i])
    
sign = { "F": -1, "Cl": -1, "Br": -1, "I": -1, "Li": 1, "Na": 1, "K": 1 }

tempdir = "vtmp"
os.makedirs(tempdir, exist_ok=True)
os.chdir(tempdir)

# Generate point charge data
minus = "minus.xvg"
with open(minus, "w") as outf:
    for i in range(start, len(xref)):
        outf.write("%10g  %10g\n" % ( xref[i], -pc[i] ))

plus = "plus.xvg"
with open(plus, "w") as outf:
    for i in range(10, len(xref)):
        outf.write("%10g  %10g\n" % ( xref[i], pc[i] ))

plots = [ "relative", "scaled", "absolute" ]
for rel in plots:
    for ion in output_data["data"].keys():
        xy = output_data["data"][ion]
        if len(xy) == 2 and (len(xy[0]) == len(xy[1])):
            pname = ("%s-%s.xvg" % ( ion, rel ))
            with open(pname, "w") as outf:
                outf.write("@ xaxis label \"Distance ($\mathrm{\AA}$)\"\n")
                relstr = ""
                if "relative" == rel:
                    relstr = " - PC"
                elif "scaled" == rel:
                    relstr = " - Scaled"
                outf.write("@ yaxis label \"ESP%s (kJ/mol e)\"\n" % relstr)
                for i in range(len(xy[0])):
                    if i < start:
                        continue
                    yy = xy[1][i]
                    if "relative" == rel:
                        yy -= sign[ion]*pc[i]
                    elif "scaled" == rel:
                        yy /= sign[ion]*pc[i]
                    outf.write("%10g  %10g\n" % ( xy[0][i], yy ))

    opts =  "-mk  -xframe 12 -yframe 6  -noshow -legend_x 0.53 -legend_y 1.05 -alfs 48 -tickfs 54 -lfs 38 -colors crimson royalblue lightcoral mediumseagreen"

    if "relative" == rel:
        opts += " -ymax 200"
        os.system("viewxvg -f F-%s.xvg Cl-%s.xvg Br-%s.xvg -label Fluoride Chloride Bromide -pdf ../anion-esp-%s.pdf %s" %
                  (rel, rel, rel, rel, opts))
        os.system("viewxvg -f Li-%s.xvg Na-%s.xvg K-%s.xvg -label Lithium Sodium Potassium -pdf ../cation-esp-%s.pdf %s" %
                  (rel, rel, rel, rel, opts))
    elif "scaled" == rel:
        os.system("viewxvg -f F-%s.xvg Cl-%s.xvg Br-%s.xvg -label Fluoride Chloride Bromide -pdf ../anion-esp-%s.pdf %s -ymax 2" %
                  (rel, rel, rel, rel, opts ))
        os.system("viewxvg -f Li-%s.xvg Na-%s.xvg K-%s.xvg -label Lithium Sodium Potassium -pdf ../cation-esp-%s.pdf %s -ymax 4" %
                  (rel, rel, rel, rel, opts ))
    else:
        os.system("viewxvg -f F-%s.xvg Cl-%s.xvg Br-%s.xvg minus.xvg -label Fluoride Chloride Bromide PC -pdf ../anion-esp-%s.pdf %s" %
                  (rel, rel, rel, rel, opts))
        os.system("viewxvg -f Li-%s.xvg Na-%s.xvg K-%s.xvg plus.xvg -label Lithium Sodium Potassium PC -pdf ../cation-esp-%s.pdf %s" %
                  (rel, rel, rel, rel, opts))
