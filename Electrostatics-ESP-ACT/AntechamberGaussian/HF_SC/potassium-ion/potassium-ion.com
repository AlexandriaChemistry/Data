%mem=12000MB
%nprocshared=16
%mem=32Gb
%chk=potassium-ion.chk
#P HF/aug-cc-pvtz Opt=tight  Pop=(MK,Hirshfeld) iop(6/33=2) iop(6/42=6)
maxdisk=128GB
            
potassium-ion

1 1
K          0.00000        0.00000        0.00000


