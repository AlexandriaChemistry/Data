%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB
            
ammonium

1 1
    H   -0.00004123    1.02353975    0.00000000
    N    0.00001806    0.00003319    0.00000000
    H   -0.96526601   -0.34134680    0.00000000
    H    0.48252818   -0.34132705   -0.83578018
    H    0.48252818   -0.34132705    0.83578018


