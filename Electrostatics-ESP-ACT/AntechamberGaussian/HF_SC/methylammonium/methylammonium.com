%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            
methylammonium

1 1
    N   -0.71078813    0.00008135   -0.06696861
    C    0.81021448   -0.00011747    0.07634089
    H   -1.02025782   -0.74795559   -0.69228706
    H   -1.04299911    0.88813974   -0.45069962
    H   -1.16442794   -0.14007795    0.83921842
    H    1.24759713    0.15104533   -0.90513070
    H    1.11672008   -0.95703408    0.48497612
    H    1.09222182    0.80615091    0.74543027


