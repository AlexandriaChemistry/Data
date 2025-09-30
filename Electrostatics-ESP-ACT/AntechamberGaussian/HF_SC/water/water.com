%mem=12000MB
%nprocshared=16
%mem=32Gb
%chk=water.chk
#P HF/aug-cc-pvtz Opt=tight  Pop=(MK,Hirshfeld) iop(6/33=2) iop(6/42=6)
maxdisk=128GB
            
water

0 1
O          0.00000        0.00000        0.11285
H         -0.00000        0.75309       -0.45140
H         -0.00000       -0.75309       -0.45140


