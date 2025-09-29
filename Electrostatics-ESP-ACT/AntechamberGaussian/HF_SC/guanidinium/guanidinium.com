%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            
guanidinium

1 1
C         -0.00730        0.02480       -0.00000
N         -1.09310       -0.75020        0.00000
N          1.20680       -0.52800       -0.00000
N         -0.13560        1.35270       -0.00000
H          2.06690        0.01150        0.00000
H          1.36460       -1.53100       -0.00000
H         -2.04060       -0.38530       -0.00000
H         -1.05590       -1.76490       -0.00000
H         -1.03290        1.82780        0.00000
H          0.65420        1.99080        0.00000


