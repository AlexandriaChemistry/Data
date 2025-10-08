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
    for i in range(start, len(xref), delta):
        f_minus.write(f"{xref[i]:10.4f}  {-pc[i]:10.6f}\n")

with open("plus.xvg", "w") as f_plus:
    labels(f_plus)
    for i in range(start, len(xref), delta):
        f_plus.write(f"{xref[i]:10.4f}  {pc[i]:10.6f}\n")

for ion, xy in output_data["data"].items():
    if len(xy) == 2 and len(xy[0]) == len(xy[1]):
        pname = f"{ion}-absolute.xvg"
        with open(pname, "w") as outf:
            labels(outf)
            for i in range(start, len(xy[0]), delta):
                outf.write(f"{xy[0][i]:10.4f}  {xy[1][i]:10.6f}\n")

opts = (
    "-mk o x + . -xframe 12 -yframe 12 -noshow "
    "-legend_x 0.5 -legend_y 0.8 "
    "-alfs 48 -tickfs 48 -lfs 38 "
#    "-colors crimson lightcoral royalblue mediumseagreen"
)

os.system(
    f"viewxvg -f F-absolute.xvg Cl-absolute.xvg Br-absolute.xvg minus.xvg "
    f"-label Fluoride Chloride Bromide PC "
    f"-save ../anion-esp-absolute.pdf {opts}"
)

os.system(
    f"viewxvg -f Li-absolute.xvg Na-absolute.xvg K-absolute.xvg plus.xvg "
    f"-label Lithium Sodium Potassium PC "
    f"-save ../cation-esp-absolute.pdf {opts}"
)

