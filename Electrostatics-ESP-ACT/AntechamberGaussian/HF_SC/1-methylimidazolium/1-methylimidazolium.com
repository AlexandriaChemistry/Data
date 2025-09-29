%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            
1-methylimidazolium

1 1
N         -0.71470       -1.44870       -0.00000
C          0.66150       -1.47200       -0.00000
C          1.10320       -0.19600       -0.00000
N          0.00220        0.66250       -0.00000
C         -1.10370       -0.14950       -0.00000
C          0.01730        2.10180        0.00000
H          1.21900       -2.39520       -0.00000
H          2.11380        0.18160        0.00000
H         -2.12810        0.19850        0.00000
H          1.04090        2.49170       -0.00000
H         -0.49150        2.49090        0.88820
H         -0.49150        2.49090       -0.88820
H         -1.32090       -2.25160        0.00000


