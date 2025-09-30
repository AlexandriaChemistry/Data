%mem=12000MB
%nprocshared=16
%mem=32Gb
%chk=methylammonium.chk
#P HF/aug-cc-pvtz Opt=tight  Pop=(MK,Hirshfeld) iop(6/33=2) iop(6/42=6)
maxdisk=128GB
            
methylammonium

1 1
N         -0.00000        0.00000       -0.70716
C          0.00000       -0.00000        0.79612
H         -0.00000        0.93845       -1.07406
H         -0.81272       -0.46922       -1.07406
H          0.81272       -0.46922       -1.07406
H         -0.88634        0.51173        1.13187
H          0.88634        0.51173        1.13187
H          0.00000       -1.02346        1.13187


