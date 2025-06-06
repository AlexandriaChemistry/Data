#!/usr/bin/env python3

import os

acmparm = {
    "ACM-esp-G":      { "ff": "esp-g.xml", "nparm": 48, "label": "GC", "target": "ESP" },
    "ACM-esp-GV":     { "ff": "esp-gv.xml", "nparm": 54, "label": "GC+PGV", "target": "ESP" },
    "ACM-elec-P":     { "ff": "coul-p.xml", "nparm": 32, "label": "PC", "target": "Elec" },
    "ACM-allelec-P":  { "ff": "all-p.xml", "nparm": 32, "label": "PC", "target": "Elec+Induc" },
    "ACM-elec-G":     { "ff": "coul-g.xml", "nparm": 48, "label": "GC", "target": "Elec" },
    "ACM-allelec-G":  { "ff": "all-g.xml", "nparm": 48, "label": "GC", "target": "Elec+Induc" },
    "ACM-elec-GV":    { "ff": "coul-gv.xml", "nparm": 54, "label": "GC+PGV", "target": "Elec" },
    "ACM-allelec-GV": { "ff": "all-gv.xml", "nparm": 54, "label": "GC+PGV", "target": "Elec+Induc" },
    "ACM-all-PG":     { "ff": "all-pg.xml", "nparm": 123, "label": "PC+GVS", "target": "Elec,Induc" }
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
        if model == "ACM-all-PG":
            os.system("grep -v a2dexp %s | grep -v hyper > koko.xml" % myff)
            os.system("mv koko.xml %s" % myff)
        
        texfn = ( "params-%s.tex" % model )
        os.system("alexandria merge_ff -ff %s -latex %s -info '%s' -merge '%s'" % ( myff, texfn, myinfo, params ) )
        prune_eem(texfn)
