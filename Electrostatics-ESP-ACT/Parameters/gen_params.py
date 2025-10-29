#!/usr/bin/env python3

import os

header = "header"

nq = 45
acmparm = {
    "PC+GV-esp3":    { "ff": "PC+GV-esp3.xml", "nparm": nq*3, "label": "PC+GV3x", "target": "ESP", header: "Non-polarizable ACT monomer-based models" },
    "PC+GV-esp4":    { "ff": "PC+GV-esp4.xml", "nparm": 2+nq*3, "label": "PC+GV4x", "target": "ESP" },
    "PC+SV-esp3":    { "ff": "PC+SV-esp3.xml", "nparm": nq*3, "label": "PC+SV3x", "target": "ESP" },
    "PC+SV-esp4":    { "ff": "PC+SV-esp4.xml", "nparm": 2+nq*3, "label": "PC+SV4x", "target": "ESP" },
    "PC-elec":       { "ff": "PC-elec.xml", "nparm": 66, "label": "PC", "target": "Elec", header: "Non-polarizable ACT dimer-based models" },
    "PC-allelec":    { "ff": "PC-allelec.xml", "nparm": 66, "label": "PC", "target": "Elec+Induc" },
    "GC-elec":       { "ff": "GC-elec.xml", "nparm": 85, "label": "GC", "target": "Elec" },
    "GC-allelec":    { "ff": "GC-allelec.xml", "nparm": 85, "label": "GC", "target": "Elec+Induc" },
    "SC-elec":       { "ff": "SC-elec.xml", "nparm": 85, "label": "SC", "target": "Elec" },
    "SC-allelec":    { "ff": "SC-allelec.xml", "nparm": 85, "label": "SC", "target": "Elec+Induc" },
    "PC+GV-elec":    { "ff": "PC+GV-elec.xml", "nparm": 106, "label": "PC+GV4", "target": "Elec" },
    "PC+GV-allelec": { "ff": "PC+GV-allelec.xml", "nparm": 106, "label": "PC+GV4", "target": "Elec+Induc" },
    "PC+SV-elec":    { "ff": "PC+SV-elec.xml", "nparm": 106, "label": "PC+SV4", "target": "Elec" },
    "PC+SV-allelec": { "ff": "PC+SV-allelec.xml", "nparm": 106, "label": "PC+SV4", "target": "Elec+Induc" },
    "PC+GS-elec":    { "ff": "PC+GS-elec.xml", "nparm": 156, "label": "PC+GS4", "target": "Elec,Induc", header: "Polarizable ACT dimer-based models" },
    "PC+GS-allelec": { "ff": "PC+GS-allelec.xml", "nparm": 156, "label": "PC+GS4", "target": "Elec+Induc" }
}

def prune_eem(tab:str):
    temp = "koko.tex"
    os.system("cp %s %s" % ( tab, temp ))
    with open(temp, "r") as inf:
        with open(tab, "w") as outf:
            foundchi = False
            for line in inf:
                skipLine = False
                if "chi" in line and not "delta" in line:
                    foundchi = True
                if foundchi:
                    if "longtable" in line:
                        foundchi = False
                    else:
                        words = line.split("&")
                        if len(words) == 3:
                            four = "4.000"
                            if four in words[1] and four in words[2]:
                                skipLine = True
                if not skipLine:
                    outf.write(line)
    
if __name__ == "__main__":
    atoms = ""
    for aa in [ "Ar", "he", "ar", "kr", "xe", "ne", "Be2+", "Ca2+", "Cs+", "He", "I-", "Kr", "Mg2+", "Ne", "Rb+", "Xe", "b", "br", "cl", "c1", "f", "h3", "h4", "h5", "h6",
                "ha", "hbr", "hcl", "hf", "hi", "ho", "hp", "hs", "i", "n1", "n2", "n3", "o1", "o3", "p1", "p2", "p3", "s1", "s2", "s3",
                "v1Cs+", "v1I-", "v1Rb+", "v2br", "v2cl", "v2core", "v2f", "v2i", "v2o1", "v3Ow" ]:
        atoms += (" %s %s_s %s_z %s_s_z " % ( aa, aa, aa, aa ) )
    params = "zeta alpha chi eta delta_eta delta_chi vs3sa a1dexp bdexp"
    for model in acmparm.keys():
        myinfo = ("Model: %s, training target: %s." % ( acmparm[model]["label"], acmparm[model]["target"] ) )
        myff   = acmparm[model]["ff"]
        os.system("alexandria edit_ff -ff ../AlexandriaFF/%s -o %s -del -a '%s'" % ( myff, myff, atoms ) ) 
        if model.startswith("PC+GS"):
            os.system("grep -v a2dexp %s | grep -v hyper > koko.xml" % myff)
            os.system("mv koko.xml %s" % myff)
        
        texfn = ( "params-%s.tex" % model )
        os.system("alexandria merge_ff -ff %s -latex %s -info '%s' -merge '%s'" % ( myff, texfn, myinfo, params ) )
        prune_eem(texfn)
