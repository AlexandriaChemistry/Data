%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            
acetate

-1 1
    C    0.00276774    0.00016340   -1.44443037
    C    0.00494962    0.00066195    0.11409143
    O    1.15398988    0.00033222    0.68194663
    O   -1.15650435   -0.00077807    0.66007789
    H    1.02001481    0.01687288   -1.83803420
    H   -0.52038137   -0.88985500   -1.80971488
    H   -0.55161635    0.87023073   -1.81103591


