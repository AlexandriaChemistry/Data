#!/usr/bin/env python3

import os
from run_gaussian import get_mols

def extract_columns(input_file, skip_lines=7)->list:
    atomq = []
    with open(input_file, 'r') as infile:
        for _ in range(skip_lines):
            next(infile)
        
        for line in infile:
            if line.strip() == 'LOOP':
                break
            
            parts = line.split()
            if len(parts) == 11 and parts[1] != 'DUMM':
                atomq.append( { "name": parts[1], "q": parts[-1]  })
    return atomq

if __name__ == "__main__":
    mols = get_mols()
    for mol in mols.keys():
        mdir = f"HF_SC/{mol}"
        if os.path.isdir(mdir):
            txml = f"{mdir}/temp.xml"
            os.system(f"gauss2molprop -n {mol} -i {mdir}/{mol}.log -fi {pdb} -o {txml} -basis aug-cc-pvtz")
            
            atomq = {}
            for method in [ "bcc", "resp" ]:
                atomq[method] = extract_columns(f"prepi/{mol}_{method}.prepi")
            
            
    os.system("alexandria edit_mp -mp */*.xml -o hf-aug-cc-pvtz.xml")
