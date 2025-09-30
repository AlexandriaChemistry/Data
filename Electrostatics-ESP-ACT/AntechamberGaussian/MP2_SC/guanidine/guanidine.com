%mem=12000MB
%nprocshared=16
%mem=32Gb
%chk=guanidine.chk
#P MP2/aug-cc-pvtz Opt=tight  density=MP2 Pop=(MK,Hirshfeld) iop(6/33=2) iop(6/42=6)
maxdisk=128GB
            
guanidine

0 1
C         -0.01906        0.12240        0.00027
N         -0.94181       -0.90385        0.07326
N          1.27860       -0.33497       -0.07815
N         -0.26270        1.35647        0.00596
H          1.94617        0.39223        0.03425
H          1.47925       -1.14265        0.46501
H         -1.88677       -0.63140       -0.05967
H         -0.70078       -1.72806       -0.42685
H         -1.24212        1.55200       -0.02185


