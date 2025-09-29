%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            
formate

-1 1
    O   -1.15890697   -0.18124493    0.00000000
    O    1.15897748   -0.18122310    0.00000000
    C   -0.00009064    0.35780328    0.00000000
    H   -0.00003980    1.49232837    0.00000000


