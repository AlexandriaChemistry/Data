# Accurate-Electrostatics-for-Biomolecular-Systems
Scripts, force field files for widely-used and Alexandria force fields (only electrostatic part), and electrostatic potential (ESP) fitting for reproducing results in the scientific article
[_Accurate Electrostatics for Biomolecular Systems through Machine Learning_].
Here's the [preprint] (https://chemrxiv.org/engage/chemrxiv/article-details/67f82e8e81d2151a02e49d35)

## Directory Layout

- `AnalyticalFitting/` - Analytical fitting for the quantum chemical ESP and electrostatic calculations for different charge models, e.g. a (positive) point charge with either one Gaussian or 1S Slater distributed charge, or a
   point charge with two Gaussian charges or a point charge with a 1S and a 2S Slater charge.
- `AnalyticalFitting/SAPT_Alkali_Halides/` - SAPT0/aug-cc-pVTZ calculations.
- `AntechamberGaussian/` - BCC and RESP calculations from Gaussian and Antechamber for ions, water, and side-chains.
- `Charge_Models/` - Reproduction of Table 1, Comparison of electrostatic energies from Alexandria and other popular charge models, e.g. ESP. Note that according to the definition of SAPT, we 
   need the sum of electrostatic and induction energies from SAPT to obtain the correct electrostatic energy. 
- `Partial-Charges/` - Reproduction of partial charges from log files- Tables S8-S13.  
- `SAPTDatabase/` - SAPT2+(CCD)-δMP2 database for training and tet.
- `Scripts4Figs/` - Reproduction of the figure 2B,S2 presented in the article.
- `AlexandriaFF/` - Trained force field files on SAPT2+(CCD)-δMP2 using Alexandria.esp-paper-gaussian.xml FF contains charges for Mulliken, Hirshfeld, CM5, and ESP.
- `ForceFields/` - Available force field files such as TIP3P, TIP4P, TIP4P-2005.xml, TIP4PEW-JC.xml, Walz2018a.xml, GAFF, and SWM4-NDP.xml. 
- `Selection/` - Train and Test sets and a script to generate random test and train compounds. 
- `Scripts4Tabs/` - Reproduction of tables presented in the article.
- `Training_RMSD_Compound/` - RMSD values with respect to SAPT for all studied compounds, including our ACT model and several widely used models. 
- `Training_RMSD_Compound/Elec/` - Training files on electorsrtatic energy from SAPT2+(CCD)-δMP for PC, GC, PC+GV, and PC+GVS models.  
- `Training_RMSD_Compound/AllElec/` - Training files on electorsrtatic and induction energies from SAPT2+(CCD)-δMP2 for PC, GC and, PC+G model.
- `Training_RMSD_Compound/PARAMETERS/` -  Parameters for the PC, GC, GC+PGV, and PC+GVS  models trained on the electrostatic energy. For details please see the ACT manual.
- `ESP-RMSD-Histogram/` - ESP at the HF/6-31G** level of theory for 5100 compounds, Fig S1.

## Requirements

- To run the scripts, you'll need Python installed on your computer, if you want to run on your own computer,
install python using e.g. [Miniconda](https://conda.io/miniconda.html) or [Anaconda](https://docs.conda.io)).

- To be able to run trainings, you need to install Alexandria Chemistry Toolkit (ACT) and its dependencies (ACT/ACT.yml).

