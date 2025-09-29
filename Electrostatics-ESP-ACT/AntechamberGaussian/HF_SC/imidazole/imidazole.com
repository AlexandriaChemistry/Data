%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            
imidazole

0 1
    C    1.14525960    0.32033900    0.00000000
    H    2.14135833    0.72109623    0.00000000
    C    0.66737273   -0.96716053    0.00000000
    H    1.21655121   -1.88930074    0.00000000
    N   -0.72964267   -0.96636240    0.00000000
    C   -1.08858766    0.31217954    0.00000000
    H   -2.09700378    0.68100836    0.00000000
    N    0.01727701    1.13430633    0.00000000
    H    0.01587641    2.13824912    0.00000000


