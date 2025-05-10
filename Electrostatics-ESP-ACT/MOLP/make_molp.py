#!/usr/bin/env python 
import subprocess


data_homo = [
    "water#water", "ammonium#chloride", "ammonium#fluoride", "ammonium#formate", "ammonium#potassium",
    "ethylammonium#formate", "methylammonium#formate", "methylammonium#acetate", "ethylammonium#ethylammonium",
    "sodium#sodium", "potassium#potassium", "chloride#chloride", "bromide#bromide", "ammonium#acetate",
    "ammonium#bromide", "ammonium#water", "ammonium#lithium", "ammonium#sodium", "ethylammonium#bromide",
    "ethylammonium#chloride", "chloride#fluoride", "formate#chloride", "bromide#fluoride", "bromide#chloride",
    "ethylammonium#methylammonium", "sodium#potassium", "acetate#acetate", "acetate#bromide", "lithium#lithium",
    "fluoride#fluoride", "lithium#potassium", "ammonium#ethylammonium", "formate#acetate", "methylammonium#methylammonium",
    "methylammonium#lithium", "methylammonium#sodium", "ammonium#ammonium", "formate#formate", "lithium#sodium",
    "ammonium#methylammonium"
]

ions = [
    "water#bromide", "water#lithium", "water#potassium", "potassium#chloride", "acetate#lithium", "acetate#sodium",
    "acetate#potassium", "acetate#water", "acetate#fluoride", "acetate#chloride", "methylammonium#fluoride",
    "methylammonium#chloride", "methylammonium#bromide", "methylammonium#water", "ethylammonium#water",
    "formate#potassium", "formate#water", "sodium#chloride", "potassium#bromide", "ethylammonium#potassium",
    "sodium#fluoride", "lithium#bromide", "water#chloride", "water#fluoride", "water#sodium", "formate#sodium",
    "formate#lithium", "ethylammonium#fluoride", "lithium#fluoride", "ethylammonium#lithium", "ethylammonium#sodium",
    "formate#fluoride", "formate#bromide", "methylammonium#potassium", "potassium#fluoride", "sodium#bromide",
    "lithium#chloride"
]

with open("homo.dat", "w") as f:
    f.write("\n".join(data_homo))

with open("het.dat", "w") as f:
    f.write("\n".join(ions))


subprocess.run([
    "./write_molprop.py", "-method", "sapt2+(ccd)dmp2", "-basis", "aug-cc-pvtz",
    "-o", "sapt2-homo-0.015.xml", "-dEmax", "0.015", "-sel", "homo.dat", "-rmax", "6"
])
subprocess.run([
    "./write_molprop.py", "-method", "sapt2+(ccd)dmp2", "-basis", "aug-cc-pvtz",
    "-o", "sapt2-het-0.08.xml", "-dEmax", "0.08", "-sel", "het.dat", "-rmax", "6"
])

with open("sapt2-homo-0.015.xml", "r") as f1, open("sapt2-het-0.08.xml", "r") as f2:
    lines1 = f1.readlines()
    lines2 = f2.readlines()


remove_block = [
    "</molecules>\n",
    '<?xml version="1.0" encoding="ISO-8859-1"?>\n',
    '<!DOCTYPE molprops.dtd PUBLIC "molprops.dtd" "molprops.dtd">\n',
    "<molecules>\n"
]

for i in range(len(lines2) - len(remove_block) + 1):
    if lines2[i:i+len(remove_block)] == remove_block:
        lines2 = lines2[:i] + lines2[i+len(remove_block):]
        break  

with open("merged.xml", "w") as f:
    f.writelines(lines1 + lines2)
