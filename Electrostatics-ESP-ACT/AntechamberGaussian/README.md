# Scripts for monomer calculations
This directory contains Hartree-Fock level calculations using the aug-cc-pVTZ basis set for side chain analogs, performed with Gaussian (located in the HF_SC folder). It also includes the command lines used to obtain BCC and RESP charges via Antechamber, using the Amber/22-nsc1-gcccuda-11.4-9.3.0-bare module, which are provided in commands.txt. The read_prepi.py script reads the charges from the generated prepi files.

The run_mbis script does MBIS calculations, it needs a development version of Psi4 to do so, since it needs to extract core and valence charges.
