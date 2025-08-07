The following XML files represent different force field models:

- `all-g.xml`: Charge model with Gaussian-distributed charges (GCs).

- `all-gv.xml`: Non-polarizable model (GC + PGV); a virtual site was added to the halide ions and potassium, each carrying a Gaussian-distributed charge.

- `all-p.xml`: Point charge (PC) model.

These force fields were trained on the sum of electrostatic and induction energy components from SAPT.

In addition:
- `all-pg-submit.xml`: A polarizable Gaussian-distributed shell (Drude particle) was added to generate the PC + GVS model were trained on the electrostatic energy from SAPT.

- `coul-g.xml`, `coul-gv.xml`, and  `coul-p.xml` were trained only on the electrostatic energy from SAPT.

- `esp-g.xml` and `esp-gv.xml` were trained only on the electrostatic potentials (ESP), where the quantum chemical ESP of monomeric compounds is employed.

- `esp-paper-gaussian.xml` includes available charge models, and `sapt-esp.xml` serve as reference XML files.

- `hf-aug-cc-pvtz.xml` serves as reference optimized strucrure for charge calculations. 
