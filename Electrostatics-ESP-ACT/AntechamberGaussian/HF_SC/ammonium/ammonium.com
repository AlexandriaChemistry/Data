%mem=12000MB
%nprocshared=16
%mem=32Gb
%chk=ammonium.chk
#P HF/aug-cc-pvtz Opt=tight  Pop=(MK,Hirshfeld) iop(6/33=2) iop(6/42=6)
maxdisk=128GB
            
ammonium

1 1
H         -0.58309        0.58309        0.58309
N          0.00000        0.00000        0.00000
H          0.58309       -0.58309        0.58309
H         -0.58309       -0.58309       -0.58309
H          0.58309        0.58309       -0.58309


