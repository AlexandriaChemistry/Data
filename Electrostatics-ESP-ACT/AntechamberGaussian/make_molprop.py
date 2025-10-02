#!/usr/bin/env python3

import json, os
from run_gaussian import get_mols

debug   = False
atoms   = "atom_names"
charges = "charges"
widths  = "valence_widths"

def read_mbis(mol:str, method:str)->list:
    js = f"../MBIS/{mol}_mbis_ps.json"
    if not os.path.exists(js):
        print(f"Cannot find {js}")
        return None
    mb = []
    with open(js, "r") as inf:
        data = json.load(inf)
        if atoms in data and charges in data and widths in data:
            for i in range(len(data[atoms])):
                mb.append( { "atom": data[atoms][i],
                             "q": data[charges][i][0],
                             "width": data[widths][i][0] } )
    return mb
            
def extract_columns(mol:str, input_file, skip_lines=7)->list:
    renum = { 
        "ammonium": { "4": 2, "5": 1, "6": 3, "7": 4, "8": 5 },
        "acetate":  { "4": 1, "5": 5, "6": 6, "7": 7, "8": 2, "9": 3, "10": 4 },
        "butanoate": { "4": 3, "5": 2, "6": 4, "7": 1, "8": 6, "9": 7, "10": 5,
                       "11": 8, "12": 9, "13": 10, "14": 11, "15": 12, "16": 13 },
        "ethylammonium": {  "4": 1, "5": 3, "6": 4, "7": 5, "8": 2, "9": 6, "10": 7,
                            "11": 8, "12": 9, "13": 10, "14": 11 },
        "formate": { "4": 1, "5": 3, "6": 4, "7": 2 },
        "guanidine": { "4": 2, "5": 7, "6": 8, "7": 1, "8": 4, "9": 9, "10": 3,
                       "11": 5, "12": 6 },
        "guanidinium": { "4": 2, "5": 7, "6": 8, "7": 1, "8": 3, "9": 5, "10": 6,
                         "11": 4, "12": 9, "13": 10 },
        "imidazole": { "4": 1, "5": 2, "6": 3, "7": 4, "8": 5, "9": 6, "10": 7,
                         "11": 8, "12": 9 },
        "imidazolium": { "4": 1, "5": 6, "6": 2, "7": 7, "8": 3, "9": 8, "10": 4,
                         "11": 9, "12": 5, "13": 10 },
        "1-methylimidazolium": { "4": 3, "5": 8, "6": 2, "7": 7, "8": 1, "9": 13, "10": 5,
                                 "11": 9, "12": 4, "13": 6, "14": 10, "15": 11, "16": 12 },
        "methylammonium": { "4": 1, "5": 3, "6": 4, "7": 5, "8": 2, "9": 6, "10": 7,
                            "11": 8 },
        "propanoate": { "4": 1, "5": 5, "6": 2, "7": 3, "8": 6, "9": 7, "10": 4,
                        "11": 8, "12": 9, "13": 10 },
        "water": { "4": 1, "5": 2, "6": 3 }
    }
    if not mol in renum:
        monatom = { "bromide": { "name": "Br", "q": "-1", "index": 1 },
                    "chloride": { "name": "Cl", "q": "-1", "index": 1 },
                    "fluoride": { "name": "F", "q": "-1", "index": 1 },
                    "lithium-ion": { "name": "Li", "q": "1", "index": 1 },
                    "sodium-ion": { "name": "Na", "q": "1", "index": 1 },
                    "potassium-ion": { "name": "K", "q": "1", "index": 1 } }
        if mol in monatom:
            return [ monatom[mol] ]
        else:
            return None
    if not os.path.exists(input_file):
        return None
    # Create table with right length
    atomq = [-1]*len(renum[mol])
    with open(input_file, 'r') as infile:
        for _ in range(skip_lines):
            next(infile)
        
        for line in infile:
            if line.strip() == 'LOOP':
                break
            
            parts = line.split()
            if len(parts) == 11 and parts[1] != 'DUMM':
                index = parts[0]
                if not index in renum[mol]:
                    print(f"Cannot find index {index} in renum[{mol}]")
                    return None
                atomq[renum[mol][index]-1] = { "name": parts[1], "q": parts[-1], "index": renum[mol][index] }
    return atomq

if __name__ == "__main__":
    repl = { "InChI=1S/H4N/h1H4": "InChI=1S/H3N/h1H3/p+1",
             "InChI=1S/C2H8N/c1-2-3/h2H2,1,3H3": "InChI=1S/C2H7N/c1-2-3/h2-3H2,1H3/p+1",
             "InChI=1S/CH6N/c1-2/h1-2H3": "InChI=1S/CH5N/c1-2/h2H2,1H3/p+1",
             "InChI=1S/CH6N3/c2-1(3)4/h2-4H2": "InChI=1S/CH6N3/c2-1(3)4/h2-4H2/q+1",
             "InChI=1S/C4H8O2/c1-2-3-4(5)6/h2-3H2,1H3,(H,5,6)/p-1": "InChI=1S/C4H9O2/c1-2-3-4(5)6/h2H2,1,3H3,(H,5,6)/q-1" }
    mols = get_mols()
    for qm in [ "HF", "MP2" ]:
        myxmls = ""
        for mol in mols.keys():
            mdir = f"{qm}_SC/{mol}"
            if os.path.isdir(mdir):
                txml   = f"{mdir}/temp.xml"
                molxml = f"{mdir}/{mol}.xml"
                myxmls += " " + molxml
                os.system(f"gauss2molprop -n {mol} -i {mdir}/{mol}.log -o {txml} -basis aug-cc-pvtz")
            
                atomq = {}
                for method in [ "bcc", "resp" ]:
                    atomq[method] = extract_columns(mol, f"prepi/{mol}_{qm}_{method}.prepi")
                    if debug:
                        print(f"mol {mol} method {method} atomq {atomq[method]}")
                mbisf  = None
                mbdata = None
                if qm == "MP2":
                    mbdata = read_mbis(mol, qm)
                    if debug:
                        print(f"mol {mol} mbdata {mbdata}")
                    if mbdata:
                        mbxml  = f"{mdir}/mbis.xml"
                        mbisf  = open(mbxml, "w")
                with open(molxml, "w") as outf:
                    with open(txml) as fd:
                        iatom = 0
                        for line in fd:
                            myline = line
                            for r in repl:
                                myline = myline.replace(r, repl[r])
                            if mbisf:
                                mbisf.write(myline)
                                if line.find("qESP") >= 0:
                                    # Fake it until you make it: use the RESP charge field
                                    # to store the MBIS charge.
                                    myq = mbdata[iatom]["q"]
                                    mbisf.write(f"      <qRESP>{myq}</qRESP>\n")

                            outf.write(myline)
                            if line.find("qMulliken") >= 0:
                                for method in [ "bcc", "resp" ]:
                                    if method in atomq and atomq[method]:
                                        m   = method.upper()
                                        if debug:
                                            print(f"mol {mol} iatom {iatom} atomq[method] {atomq[method]}")
                                        myq = atomq[method][iatom]["q"]
                                        outf.write(f"      <q{m}>{myq}</q{m}>\n")
                                iatom += 1
                    os.unlink(txml)
                if mbisf:
                    mbisf.close()
        if qm == "MP2":
            os.system(f"alexandria edit_mp -mp {qm}_SC/*/mbis.xml -o mbis_ccsd.xml")
        os.system(f"alexandria edit_mp -mp {myxmls} -o {qm}-aug-cc-pvtz.xml")
