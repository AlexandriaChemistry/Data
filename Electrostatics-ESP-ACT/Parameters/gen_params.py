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

if __name__ == "__main__":
    for model in acmparm.keys():
        myinfo = ("Model: %s, training target: %s." % ( acmparm[model]["label"], acmparm[model]["target"] ) )
        os.system("alexandria merge_ff -ff ../AlexandriaFF/%s -latex params-%s -info '%s'" % ( acmparm[model]["ff"], model, myinfo ) )
