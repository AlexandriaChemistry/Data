#!/bin/sh


for c in "COULOMB-ESP.xvg" "COULOMB-MBIS-S.xvg" "COULOMB-PC+GV-elec.xvg" "COULOMB-PC+GS-elec.xvg" "COULOMB-PC-elec.xvg" "COULOMB-GC-elec.xvg" "COULOMB-Mulliken.xvg" "COULOMB-Hirshfeld.xvg" "COULOMB-BCC.xvg" "COULOMB-CM5.xvg" "COULOMB-MBIS.xvg" "COULOMB-RESP.xvg" 
do
    sed -i '' 's/label \"Alexandria\"/label \"Energy\"/g' $c
done

viewxvg -f  COULOMB-MBIS.xvg COULOMB-MBIS-S.xvg COULOMB-PC+GV-esp3.xvg COULOMB-PC+SV-esp3.xvg COULOMB-PC+GV-esp4.xvg COULOMB-PC+SV-esp4.xvg  -label  MBIS MBIS-S  PC+GV3x PC+SV3x   PC+GV4x PC+SV4x  -ls None None -mk o +  -res   -panel -sharelabel  -xframe 16 -yframe 16 -save coulomb-monomer.pdf -alfs 28 -lfs 24 -tickfs 24 -color royalblue crimson -sfx -0.15 -noshow

viewxvg -f  COULOMB-Mulliken.xvg COULOMB-Hirshfeld.xvg COULOMB-CM5.xvg COULOMB-BCC.xvg COULOMB-ESP.xvg COULOMB-RESP.xvg -label Mulliken Hirshfeld  CM5 BCC ESP RESP -ls None None -mk o +  -res   -panel -sharelabel  -xframe 16 -yframe 16 -save coulomb-legacy.pdf -alfs 28 -lfs 24 -tickfs 24 -color royalblue crimson -sfx -0.15 -noshow

viewxvg -f  COULOMB-PC-elec.xvg COULOMB-GC-elec.xvg COULOMB-SC-elec.xvg COULOMB-PC+GV-elec.xvg COULOMB-PC+SV-elec.xvg COULOMB-PC+GS-elec.xvg  -label  PC GC SC PC+GV4 PC+SV4 PC+GS4 -ls None None -mk o +  -res   -panel -sharelabel  -xframe 16 -yframe 16 -save coulomb-dimer.pdf -alfs 28 -lfs 24 -tickfs 24 -color royalblue crimson -sfx -0.15 -noshow
