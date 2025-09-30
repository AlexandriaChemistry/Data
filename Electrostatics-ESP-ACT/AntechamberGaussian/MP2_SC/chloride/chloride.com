%mem=12000MB
%nprocshared=16
%mem=32Gb
%chk=chloride.chk
#P MP2/aug-cc-pvtz Opt=tight  density=MP2 Pop=(MK,Hirshfeld) iop(6/33=2) iop(6/42=6)
maxdisk=128GB
            
chloride

-1 1
Cl         0.00000        0.00000        0.00000


