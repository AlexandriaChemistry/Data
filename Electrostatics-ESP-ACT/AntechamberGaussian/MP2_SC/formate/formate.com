%mem=12000MB
%nprocshared=16
%mem=32Gb
%chk=formate.chk
#P MP2/aug-cc-pvtz Opt=tight  density=MP2 Pop=(MK,Hirshfeld) iop(6/33=2) iop(6/42=6)
maxdisk=128GB
            
formate

-1 1
O          0.00000        1.11548       -0.20564
O         -0.00000       -1.11548       -0.20564
C          0.00000        0.00000        0.31006
H          0.00000       -0.00000        1.42988


