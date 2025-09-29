%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            
water

0 1
O          0.00000       -0.00000        0.07800
H          0.00000        0.75700       -0.50800
H          0.00000       -0.75700       -0.50800


