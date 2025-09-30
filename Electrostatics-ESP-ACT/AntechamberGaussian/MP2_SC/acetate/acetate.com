%mem=12000MB
%nprocshared=16
%mem=32Gb
%chk=acetate.chk
#P MP2/aug-cc-pvtz Opt=tight  density=MP2 Pop=(MK,Hirshfeld) iop(6/33=2) iop(6/42=6)
maxdisk=128GB
            
acetate

-1 1
C         -1.33628       -0.00002       -0.00511
C          0.20970        0.00000       -0.00957
O          0.74312       -1.11089        0.00175
O          0.74308        1.11092        0.00175
H         -1.72751       -0.89026       -0.48403
H         -1.67508        0.00016        1.02842
H         -1.72753        0.89004       -0.48434


