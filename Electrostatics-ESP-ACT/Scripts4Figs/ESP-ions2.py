#!/usr/bin/env python3

import os
import json
import numpy as np
import sys

sys.path.append('../Scripts4Tabs')

from walz2018a_SAPT0_P1S2S import one_4pi_eps0


input_json = 'output_4_100.json'
with open(input_json, 'r') as json_f:
    output_data = json.load(json_f)

start = 80
delta = 5

xref = output_data["data"]["Br"][0]

pc = [0] * start
for i in range(start, len(xref)):
    if abs(xref[i]) > 1e-8:
        pc.append(one_4pi_eps0 / xref[i])
    else:
        pc.append(0.0)

sign = {"F": -1, "Cl": -1, "Br": -1, "I": -1, "Li": 1, "Na": 1, "K": 1}

tempdir = "vtmp"
os.makedirs(tempdir, exist_ok=True)
os.chdir(tempdir)

def labels(outf):
    outf.write('@ xaxis label "Distance ($\\mathrm{\\AA}$)"\n')
    outf.write('@ yaxis label "ESP (kJ/mol e)"\n')

with open("minus.xvg", "w") as f_minus:
    labels(f_minus)
    leg = [ "F-", "Cl-", "Br-", "-PC" ]
    for i in range(len(leg)):
        f_minus.write("@ s%d legend \"%s\"\n" % ( i, leg[i] ))
    for i in range(start, len(xref), delta):
        f_minus.write("%10.4f  %10.6f  %10.6f  %10.6f  %10.6f\n" %
                      ( xref[i],
                        output_data["data"]["F"][1][i],
                        output_data["data"]["Cl"][1][i],
                        output_data["data"]["Br"][1][i],
                        -pc[i]) )

with open("plus.xvg", "w") as f_plus:
    labels(f_plus)
    leg = [ "Li+", "Na+", "K+", "PC" ]
    for i in range(len(leg)):
        f_plus.write("@ s%d legend \"%s\"\n" % ( i, leg[i] ))
    for i in range(start, len(xref), delta):
        f_plus.write("%10.4f  %10.6f  %10.6f  %10.6f  %10.6f\n" %
                      ( xref[i],
                        output_data["data"]["Li"][1][i],
                        output_data["data"]["Na"][1][i],
                        output_data["data"]["K"][1][i],
                        pc[i]) )


opts = (
    "-mk o x + . -xframe 16 -yframe 16 -noshow -panel "
    "-legend_x 0.7 -legend_y 0.99 -ls solid "
    "-alfs 48 -tickfs 48 -lfs 36 "
    "-colors crimson black royalblue mediumseagreen   crimson black royalblue mediumseagreen "
)

os.system(
    f"plotxvg -f minus.xvg plus.xvg "
    f"-save ../ion-esp-absolute.pdf {opts}"
)


