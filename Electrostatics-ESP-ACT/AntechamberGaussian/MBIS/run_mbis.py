#!/usr/bin/env python3
#SBATCH -t 48:00:00
#SBATCH -c 64
#SBATCH -p CLUSTER-AMD
import os, sys, json
import numpy as np
import psi4 as psi4
import psi4.driver.p4util as p4util

def fetch_mbis_array(wfn, name):
    try:
        return np.array(wfn.array_variable(name))
    except Exception:
        return np.array(psi4.core.variable(name))

def xyz2geometry(xyzfn:str, q:int, extra:list):
    gstr = f"{q} 1\n"
    atoms = []
    with open(xyzfn, "r") as inf:
        lines = inf.readlines()
        if len(lines) > 2:
            natom = int(lines[0].strip())
            for k in range(natom):
                gstr += lines[k+2]
                atoms.append(lines[k+2].split()[0])
    for k in range(len(extra)):
        gstr += extra[k].strip() + "\n"
    print(gstr)
    return psi4.geometry(gstr), atoms

def set_psi4_opts(compound:str, ncores:int):
    psi4.core.set_num_threads(ncores)
    psi4.set_options({'guess': 'read', 'reference': 'rhf', 'cachelevel': 1, 'print': 1, 'damping_percentage': 20, 
                      'mbis_radial_points': 199, 'mbis_spherical_points': 770, 'max_radial_moment': 4 })

    psi4.set_memory(ncores*3500000*1000)
    psi4_io = psi4.core.IOManager.shared_object()
    tmpdir = "/scratch"
    psi4_io.set_default_path(tmpdir)
    psi4.core.set_output_file(compound + ".log", False)

def store_mbis(energy:float, wfn, atoms:list, jsonfn:str):
    mbis_props = { "charges":         'MBIS_CHARGES',
                   "dipoles":         'MBIS_DIPOLES',
                   "quadrupoles":     'MBIS_QUADRUPOLES',
                   "octupoles":       'MBIS_OCTUPOLES',
#                   "valence_charges": 'MBIS_VALENCE CHARGES',
                   "valence_widths":  'MBIS_VALENCE WIDTHS',
                   #"volume_ratios":  'MBIS_VOLUME_RATIOS',
                   "radial_moments": 'MBIS_RADIAL_MOMENTS_<R^3>'
                  }
    # The command below does not work for unclear reason related to symmetry
    #p4util.free_atom_volumes(wfn)
    oeprops = psi4.core.OEProp(wfn)

    mydict = { "energy": energy, "atoms": atoms }
    for mbp in mbis_props:
        oeprops.add(mbis_props[mbp])
    oeprops.compute()
    psi4.core.print_variables()

    for mbp in mbis_props:
        mydict[mbp] = mydict[mbp] = fetch_mbis_array(wfn, mbis_props[mbp].replace("_", " ")).tolist()

    with open(jsonfn, 'w') as fout:
        json.dump(mydict, fout, indent=2)

def run_one(compound:str, charge:int, method:str, basis:str):
    geom, atoms = xyz2geometry(f"../../xyz/{compound}.xyz", charge, [])
    if basis:
        psi4.basis_helper("""
assign %s
""" % basis)
    else:
        psi4.basis_helper("""
assign aug-cc-pvtz
""")

    # We need gradient rather than energy to make sure we have the density.
    # It will be more expensive of course.
    energy, scfwfn = psi4.energy(method, molecule=geom, return_wfn=True)
    store_mbis(energy, scfwfn, atoms, "scf.json")
    grad, wfn      = psi4.gradient(method, molecule=geom, properties_origin=["NUCLEAR_CHARGE"], return_wfn=True)
    store_mbis(energy, wfn, atoms, method+".json")

if __name__ == "__main__":
    if len(sys.argv) < 4 :
        sys.exit("Usage: %s compound charge method [ basis ]" % sys.argv[0])
    ncores = 32
    SVAR = "SLURM_CPUS_PER_TASK"
    if SVAR in os.environ:
        ncores = int(os.environ[SVAR])
    set_psi4_opts(sys.argv[1], ncores)
    basis = None
    if len(sys.argv) == 5:
        basis = sys.argv[4]
    run_one(sys.argv[1], int(sys.argv[2]), sys.argv[3], basis)

