#!/bin/sh

viewxvg -f COULOMB-PC-elec.xvg COULOMB-GC-elec.xvg  COULOMB-PC+GV-elec.xvg -label PC GC PC+GV -mk + v o x -ls None -res -legend_x 0.6 -save figS13.pdf -alfs 28 -lfs 28 -tickfs 24 -noshow

