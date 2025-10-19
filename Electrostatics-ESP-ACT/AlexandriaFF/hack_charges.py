#!/usr/bin/env python3

import json, os, sys
json_dir = "../AntechamberGaussian/MBIS/"

def get_q(molname:str)->list:
    json_path = os.path.join(json_dir, f"{molname}/CCSD.json")
    if not os.path.exists(json_path):
        print(f"no JSON for {molname}, skipping...")
        return []

    with open(json_path) as f:
        jdata = json.load(f)

    charges = jdata.get("charges", [])
    return charges

with open("MP2-MBIS.xml", "w") as outf:
    with open("MP2-aug-cc-pvtz.xml", "r") as inf:
        molname = None
        charges = []
        qindex  = 0
        for line in inf:
            outf.write(line)
            if line.find("molname") >= 0:
                words = line.split("=")
                molname = words[1][1:-3]
                charges = get_q(molname)
                print("Found molecule %s with %d charges" % ( molname, len(charges) ) )
                qindex = 0
            elif line.find("qRESP") >= 0:
                if len(charges) > 0 and qindex < len(charges):
                    outf.write("        <qMBIS>%g</qMBIS>\n" % ( charges[qindex][0] ) )
                    qindex += 1
                else:
                    print("Consistency error for %s" % molname)
                    continue

