# Beyond Partitioning: Using Force Field Science to Evaluate Electrostatics Models
Scripts, force field files for widely-used and Alexandria force fields
(only electrostatic part), and electrostatic potential (ESP) fitting for reproducing results in the scientific article
[_Accurate Electrostatics for Biomolecular Systems through Machine Learning_].
Here's an oldish [preprint](https://chemrxiv.org/engage/chemrxiv/article-details/67f82e8e81d2151a02e49d35)

---

## Scope and Purpose

This repository supports:

- Reproducibility of all figures and tables reported in the associated publication  
- Comparative evaluation of electrostatics models across multiple charge models  
- Analytical and quantum-chemical fitting of electrostatic potentials  
- Training and validation of electrostatics-focused force-field parameters  

Each subdirectory contains its own README file with additional details.

---

## Citation

If you find our scripts or data useful, we would greatly appreciate it if you cite our work and consider starring the repository.

---

## Directory Structure

### Core Methodology and Data

- `AnalyticalFitting/`  
  Analytical fitting of quantum-chemical ESPs and electrostatic energies for different charge models, including point charges with Gaussian or Slater-type distributions and mixed representations (e.g. PC+Gaussian, PC+Slater).

- `AnalyticalFitting/SAPT_Alkali_Halides/`  
  SAPT0/aug-cc-pVTZ interaction energy calculations for alkali-halide systems.

- `AntechamberGaussian/`  
  BCC and RESP charge derivations using Gaussian and Antechamber for ions, water, and biomolecular side chains.

---

### Reproduction of Published Results

- `Charge_Models/`  
  Reproduction of tables and figures comparing electrostatic energies obtained with Alexandria and widely used charge models (e.g. ESP and MBIS).  
  Includes reproduction of Tables S1–S17.

- `Scripts4Figs/`  
  Scripts used to reproduce Figure 2 and Figures S2-S11 from the publication.

- `Scripts4Tabs/`  
  Scripts used to reproduce tables reported in the manuscript and  bar graphs from the main paper (Figures 4 and 5).

---

### Force Fields and Parameters

- `AlexandriaFF/`  
  Trained Alexandria force fields based on SAPT2+(CCD)-δMP2 reference data, along with MBIS and MBIS-S force fileds.  

- `ForceFields/`  
  Reference force-field files, including:
  - CharmmDrude
  - TIP3P, TIP4P, TIP4P-2005, TIP4P-EW  
  - Walz2018a  
  - GAFF  
  - SWM4-NDP  

- `Parameters/`  
  Trained parameters for PC, GC, PC+GV, and PC+GS models.  
  Models are trained on electrostatic energies or combined electrostatic and induction energies from SAPT2+(CCD)-δMP2.  
  See the ACT manual for additional details.

---

### Data Sets and Analysis

- `Selection/`  
  Training and test sets, along with scripts for generating randomized training and test compounds.

- `ESP-RMSD-Histogram/`  
  Electrostatic potentials computed at the HF/6-31G** level of theory for 5100 compounds (Figure 3).

---

## Requirements

- To run the scripts, you'll need Python installed on your computer, if you want to run on your own computer,
install python using e.g. [Miniconda](https://conda.io/miniconda.html) or [Anaconda](https://docs.conda.io)).

- To be able to run analyses, you need to install [Alexandria Chemistry Toolkit](https://github.com/AlexandriaChemistry/ACT).

