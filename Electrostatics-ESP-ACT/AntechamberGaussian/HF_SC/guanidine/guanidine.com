%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            
guanidine

0 1
C          0.00000        0.12040        0.00000
N         -1.08800       -0.72790        0.00000
N          1.21280       -0.52220        0.00000
N         -0.03250        1.40140        0.00000
H          2.03610        0.04880        0.00000
H          1.30490       -1.51800        0.00000
H         -2.01240       -0.34510        0.00000
H         -0.99370       -1.72580        0.00000
H         -0.98090        1.75820        0.00000


