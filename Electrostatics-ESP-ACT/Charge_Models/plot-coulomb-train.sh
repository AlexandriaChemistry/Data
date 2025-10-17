#!/bin/sh

viewxvg -f  COULOMB-ESP.xvg COULOMB-MBIS-S.xvg COULOMB-PC-elec.xvg COULOMB-GC-elec.xvg COULOMB-PC+GV-elec.xvg COULOMB-PC+GS-elec.xvg -label  ESP MBIS-S  PC-elec GC-elec PC+GV-elec PC+GS-elec -ls None -mk o + x -res   -panel -sharelabel  -xframe 16 -yframe 16 -save coulomb-train.pdf -alfs 28 -lfs 24 -tickfs 24 -color royalblue crimson -sfx -0.15

