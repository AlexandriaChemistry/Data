#!/usr/bin/env python

import glob, os, sys

compounds_of_interest = [
    "ethylammonium",
    "ammonium",
    "methylammonium",
    "acetate",
    "formate",
    "water"
]

def extract_data_from_log(log_files):
    data = {
        "ESP":            {compound: [] for compound in compounds_of_interest},
        "ACM-elec-P":     {compound: [] for compound in compounds_of_interest},
        "ACM-allelec-P":  {compound: [] for compound in compounds_of_interest},
        "ACM-elec-G":     {compound: [] for compound in compounds_of_interest},
        "ACM-allelec-G":  {compound: [] for compound in compounds_of_interest},
        "ACM-elec-GV":    {compound: [] for compound in compounds_of_interest},
        "ACM-allelec-GV": {compound: [] for compound in compounds_of_interest},
        "ACM-all-PG":     {compound: [] for compound in compounds_of_interest}
    }

    for log_file in log_files:
        if not os.path.isfile(log_file):
            print(f"File does not exist: {log_file}")
            continue

        file_type = log_file[:-4]
        if not file_type in data:
            print("Ignoring unknown log file %s" % log_file)
            continue

        print(f"Processing {file_type} file: {log_file}")

        try:
            with open(log_file, 'r') as file:
                read_data = False
                current_compound = None
                core_shell_dict = {}

                for line in file:
                    line = line.strip()

                    if "Name:" in line:
                        for compound in compounds_of_interest:
                            c2 = ( "%s#%s" % ( compound, compound ) )
                            if f"Name: {c2}" in line:
                                current_compound = compound
                                read_data = True
                                break
                        if current_compound:
                            continue  

                    if "EPOT" in line:
                        read_data = False
                    elif read_data:
                        columns = line.split()
                        if len(columns) > 4 and columns[0].isdigit():
                            aindex     = 1
                            atom_type  = columns[2] + ( "-%d" % aindex )
                            acm_value  = columns[3]
                            for mypart in  data[file_type][current_compound]:
                                if mypart["type"] == atom_type:
                                    aindex   += 1
                                    atom_type = columns[2] + ( "-%d" % aindex )
                            data[file_type][current_compound].append( { "type": atom_type, "value": acm_value } )
                            
        except Exception as e:
            print(f"Error processing file {log_file}: {e}")

    return data

def save_data_as_latex(data):
    combined_output_file = 'particle_charges.tex'

    with open(combined_output_file, 'w') as file:
        for compound in compounds_of_interest:
            file.write("\\begin{sidewaystable}\n")
            file.write("\\centering\n")
            file.write(r"\caption{Partial charges for " + compound + " from ESP and from ACT models, point charge (PC), Gaussian charge (GC), point core+Gaussian vsite (GC+PGV), and point charge + Gaussian vsite and shell (PC+GVS).  Partial charges for the PC, GC, and GC+PGV models trained on either electrostatic energy (e) or the sum of the electrostatic and induction energy (ei) from the SAPT2+(CCD)-$\\delta$MP2 method with the aug-cc-pVTZ basis set are reported. Partial charges for the PC+GVS model, trained on the electrostatic and induction energies are also provided.}")
            file.write("\n")
#            file.write("\\hspace{-1cm}\n")
            file.write("\\begin{tabular}{lcccccccc}\n")
            file.write("\\hline\n")
            file.write(fr"Particle & ESP & PC$_{{e}}$ & PC$_{{ei}}$ & GC$_{{e}}$ & GC$_{{ei}}$ & GC+PGV$_{{e}}$ & GC+PGV$_{{ei}}$ & PC+GVS \\\\")
            file.write("\n")
            file.write("\\hline\n")

            # Loop over particles, best to use the polarizable model as reference
            npart = len(data["ACM-all-PG"][compound])
            total = {}
            for method in data.keys():
                total[method] = 0
            # We use homodimers for reference and therefore only one half of the atoms is sufficient.
            for ipart in range(int(npart/2)):
                particle = data["ACM-all-PG"][compound][ipart]
                mystr = particle["type"].replace("_", "\\_")
                for method in data.keys():
                    value = None
                    if not compound in data[method]:
                        sys.exit("Could not find compound %s for method %s" % ( compound, method ))
                    for mypart in data[method][compound]:
                        if mypart["type"] == particle["type"]:
                            try:
                                fval = float(mypart["value"])
                                total[method] += fval
                                value = str(round(fval, 4))
                            except ValueError:
                                sys.exit("Don't know how to treat charge %s for type %s" % ( mypart["value"], mypart["type"] ))
                    mystr += " & "
                    if value:
                        mystr += value
                file.write("%s\\\\\n" % mystr)
                
            file.write("\\hline\n")
            file.write("Total")
            for method in data.keys():
                file.write(" & %g" % round(total[method],2) )
            file.write("\\\\\n")
            file.write("\\hline\n")
            file.write("\\end{tabular}\n")
            file.write("\\end{sidewaystable}\n")

    print("LaTeX file generated successfully!")

log_files = glob.glob("ACM*.log")
log_files.append("ESP.log")

data = extract_data_from_log(log_files)
save_data_as_latex(data)
