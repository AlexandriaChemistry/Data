#!/usr/bin/env python

import os, subprocess

debug = True

xml_files = []
tmpdir = "xml.tmp"
os.makedirs(tmpdir, exist_ok=True)
with open("distance.csv", "r") as inf:
    for line in inf:
        if not line.startswith("#"):
            words = line.strip().split(",")
            if len(words) == 2:
                dimer = words[0]
                mindist = words[1]
                if debug:
                    print("dimer %s mindist %s" % ( dimer, mindist ))
                xml_file = tmpdir + "/" + dimer + ".xml"
                dat_file = tmpdir + "/" + dimer + ".dat"
                with open(dat_file, "w") as outf:
                    outf.write("%s\n" % dimer)
                subprocess.run(["./write_molprop.py", "-method", "sapt2+3(ccd)dmp2", "-basis", "aug-cc-pvtz", "-o", xml_file, "-dEmax", "1", "-sel", dat_file, "-rmax", "7", "-rmin", mindist], check=True)
                if os.path.exists(xml_file):
                    xml_files.append(xml_file)


if len(xml_files) > 0:
    subprocess.run(["alexandria", "edit_mp", "-mp", *xml_files, "-o", "sapt2+3-tz.xml"], check=True)

