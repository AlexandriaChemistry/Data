#!/usr/bin/env python

import glob, os, sys, json, xmltodict

debug    = False
BOHR     = 0.0529177
json_dir = "../AntechamberGaussian/MBIS"

compounds_of_interest = [
    "ammonium",
    "methylammonium",
    "ethylammonium",
    "formate",
    "acetate",
    "propanoate",
    "butanoate",
    "guanidinium",
    "imidazolium",
    "water",
    "lithium-ion",
    "sodium-ion",
    "potassium-ion",
    "fluoride",
    "chloride",
    "bromide"
    ]

log_files = {
    "ESP":        { "label": "ESP", "ncol": 1, "ff": "coul-p" },
    "PC-elec":    { "label": "PC$_{e}$", "ncol": 1 },
    "PC-allelec": { "label": "PC$_{ei}$", "ncol": 1 },
    "MBIS-S":     { "label": "MBIS-S", "ncol": 2, "ff": "P+S_updated", "mp": "MP2-aug-cc-pvtz_updated" },
    "PC+SV-elec": { "label": "PC+SV4$_{e}$", "ncol": 2 },
    "PC+GV-elec": { "label": "PC+GV4$_{e}$", "ncol": 2 },
    "PC+GS-elec": { "label": "PC+GS4$_{e,i}$", "ncol": 2 },
}

def generate_logs():
    qfn = f"../AlexandriaFF/MP2-aug-cc-pvtz"
    for lf in log_files:
        print(f"Running command for {lf}")
        myff = lf
        if "ff" in log_files[lf]:
            myff = log_files[lf]["ff"]
        mymp = "../AlexandriaFF/MP2-aug-cc-pvtz.xml"
        if "mp" in log_files[lf]:
            mymp = f"../AlexandriaFF/{log_files[lf]['mp']}.xml"

        base_command = f"alexandria train_ff -nooptimize -g {lf} -sel ../Selection/monomer.dat -mp {mymp} -ff ../AlexandriaFF/{myff} "

        if "elec" in lf:
            mycmd = base_command + f" -charges {qfn} -qalg SQE"
        elif "esp" in lf:
            mycmd = base_command + f" -charges {qfn}_updated -qalg ESP"
        elif lf == "MBIS":
            mycmd = base_command + f" -qqm qMBIS -charges ../AlexandriaFF/MP2-MBIS.xml -qalg Read "
        elif lf == "MBIS-S":
            mycmd = base_command + " -qqm None "
        else:
            mycmd = base_command + f" -qqm q{lf} -charges {qfn} -qalg Read "

        
        os.system(mycmd)

def simplify(atype:str)->str:
    newrep = { "hn": [ "hna1", "hna2", "hne1", "hne2", "hne3",
                       "hnm1", "hnm2", "hnm3", "hni", "hng" ],
               "hn_s": [ "v1hn", "v1hna1", "v1hna2",
                         "v1hne1", "v1hne2", "v1hne3",
                         "v1hnm1", "v1hnm2", "v1hnm3", "v1hni", "v1hng" ],
               "h1": [  "h1e1", "h1e2", "h1e3",
                        "h1m1", "h1m2", "h1m3" ],
               "h1_s": [ "v1h1", "v1h1e1", "v1h1e2", "v1h1e3",
                         "v1h1m1", "v1h1m2", "v1h1m3" ],
               "h2": [ "h2f" ],
               "h2_s": [ "v1h2", "v1h2f" ],
               "h4": [ "h4i" ],
               "h4_s": [ "v1h4i", "v1h4" ],
               "h5": [ "h5i" ],
               "h5_s": [ "v1h5i", "v1h5" ],
               "hc": [ "hce1", "hce2", "hce3",
                       "hca1", "hca2", "hca3" ],
               "hc_s": [ "v1hca1", "v1hca2", "v1hca3", "v1hc" ],
               "hw_s": [ "v1hw" ],
               "c2": [ "c2a", "c2f", "c2i1", "c2i2", "c2g" ],
               "c2_s": [ "v1c2f", "v1c2a", "v1c2i1", "v1c2i2", "v1c2g", "v1c2" ],
               "c3": [  "c3a", "c3e1", "c3e2", "c3m" ],
               "c3_s": [ "v1c3a", "v1c3", "v1c3m", "v1c3e1", "v1c3e2"  ],
               "n2": [ "n2i", "n2g" ],
               "n2_s": [ "v1n2i", "v1n2g", "v1n2" ],
               "n4": [ "n4a", "n4e", "n4m" ],
               "n4_s": [  "v1n4", "v1n4a", "v1n4e", "v1n4m" ],
               "o2": [ "o2f", "o2a1", "o2a2" ],
               "o2_s": [ "v1o2", "v1o2f", "v1o2a1", "v1o2a2" ],
               "ow_s": [ "v1ow" ],
               "Li+_s": [ "v1Li+" ],
               "Na+_s": [ "v1Na+" ],
               "K+_s": [ "v1K+" ],
               "F-_s": [ "v1F-" ],
               "Cl-_s": [ "v1Cl-" ],
               "Br-_s": [ "v1Br-" ]
              }
    for key in newrep:
        if atype in newrep[key]:
            return key
    return atype

def get_zeta_mbiss(compound:str):
    zeta = []
    fn = f"{json_dir}/{compound}/CCSD.json"
    if os.path.exists(fn):
        with open(fn, "r") as inf:
            data = json.load(inf)
        atoms = "atoms"
        zzz   = "valence_widths"
        if atoms in data and zzz in data:
            for a in range(len(data[atoms])):
                name = data[atoms][a]
                zeta.append({ "type": name, "zeta": 1/(2*BOHR*float(data[zzz][a][0])) })
                zeta.append({ "type": name + "_s", "zeta": 0 })
    return zeta

def get_zeta(model:str, compound:str):
    zeta = {}
    bcc  = []
    eem  = []
    if log_files[model]["ncol"] == 2:
        fn = f"../AlexandriaFF/{model}.xml"
        if os.path.exists(fn):
            with open(fn, "r") as inf:
                doc = xmltodict.parse(inf.read())
                ac = "alexandria_chemistry"
                dac = doc[ac]
                inter = "interaction"
                plist = 'parameterlist'
                ident = "@identifier"
                param = "parameter"
                val   = '@value'
                                        
                for ii in doc[ac][inter]:
                    mytype = ii["@type"]
                    if mytype in [ "COULOMB", "BONDCORRECTIONS", "ELECTRONEGATIVITYEQUALIZATION" ]:
                        for pp in ii:
                            if pp == plist:
                                for entry in ii[pp]:
                                    if ident in entry:
                                        for param in entry:
                                            if val in entry[param]:
                                                ptype = entry[param]["@type"]
                                                value = entry[param][val]
                                                if mytype == "COULOMB":
                                                    idl = entry[ident][:-2]
                                                    if not model == "MBIS-S":
                                                        idl = simplify(idl)
                                                    zeta[idl] = float(value)
                                                elif mytype == "BONDCORRECTIONS":
                                                    if not entry[ident] in bcc:
                                                        bcc[entry[ident]] = {}
                                                    bcc[entry[ident]][ptype] = value
                                                else:
                                                    if not entry[ident] in eem:
                                                        eem[entry[ident]] = {}
                                                    eem[entry[ident]][ptype] = value
                        
    return zeta#, eem, bcc
        
def extract_data_from_log():
    data = {}
    for file_type in log_files:
        log_file        = f"{file_type}.log"
        data[file_type] = {compound: [] for compound in compounds_of_interest}
        if not os.path.isfile(log_file):
            print(f"File does not exist: {log_file}")
            continue

        try:
            with open(log_file, 'r') as file:
                read_data = False
                current_compound = None
                core_shell_dict = {}

                for line in file:
                    line = line.strip()

                    if "Name:" in line:
                        for compound in compounds_of_interest:
                            if f"Name: {compound}" in line:
                                current_compound = compound
                                read_data = True
                                print("Will read charges for %s from %s" % ( compound, log_file ))
                                break
                        if current_compound:
                            continue  

                    if "MP2_SC" in line:
                        read_data = False
                    elif read_data:
                        columns = line.split()
                        if len(columns) > 4 and columns[0].isdigit():
                            aindex     = 1
                            atom_type  = columns[2]
                            if not file_type == "MBIS-S":
                                atom_type = simplify(atom_type)
                            acm_value  = columns[3]
                            for mypart in data[file_type][current_compound]:
                                if mypart["type"] == atom_type:
                                    aindex   += 1
                            data[file_type][current_compound].append( { "type": atom_type, "aindex": aindex, "value": acm_value } )
            if debug:
                print(f"file {file_type} data {data[file_type]['guanidinium']}")
        except Exception as e:
            print(f"Error processing file {log_file}: {e}")

    return data

def save_data_as_latex(data):
    combined_output_file = 'particle_charges.tex'

    zeta = {}
    for method in data.keys():
        # Do we have zeta's for this method?
        if log_files[method]["ncol"] == 2:
            zeta[method] = {}
            for compound in compounds_of_interest:
                if method == "MBIS-S":
                    zeta[method][compound] = get_zeta_mbiss(compound)
                else:
                    zeta[method][compound] = get_zeta(method, compound)

    with open(combined_output_file, 'w') as file:
        for compound in compounds_of_interest:
            file.write("\\begin{sidewaystable}\n")
            file.write("\\centering\n")
            file.write(r"\caption{Partial charges q (e) and screening widths $\zeta$ (1/nm) for " + compound + " from ESP, MBIS-S and ACT models. First line, atom, second line (X) corresponds to shell (PC+GS4) or virtual site (other models).}")
            file.write("\n")
            file.write("\\label{tab:q%s}\n" % compound)
            file.write("\\begin{tabular}{l")
            for c in data.keys():
                file.write("c")
                # If we print zeta we need to add another column
                if log_files[c]["ncol"] == 2:
                    file.write("c")
            file.write("}\n")
            file.write("\\hline\n")
            file.write("Particle ")
            for c in data:
                file.write(" & \\multicolumn{%d}{c}{%s}" % ( 
                    log_files[c]["ncol"], log_files[c]["label"] ) )
            file.write("\\\\\n")
            for c in data:
                if log_files[c]["ncol"] == 1:
                    file.write(" & q")
                else:
                    file.write(" & q & $\\zeta$" )
            file.write("\\\\\n")
            file.write("\\hline\n")

            # Loop over particles, best to use the polarizable model as reference
            PCGS = "PC+GS-elec"
            npart = len(data[PCGS][compound])
            total = {}
            for method in data.keys():
                total[method] = 0
            mystrings = []
            for ipart in range(npart):
                # Assume PC+GC has all particles
                particle = data[PCGS][compound][ipart]
                ptype    = particle["type"]
                if ptype.find("_") > 0:
                    mystr = [ ptype.split("_")[0] + " X" ]
                else:
                    mystr    = [ ptype ]
                for method in data.keys():
                    if not compound in data[method]:
                        sys.exit("Could not find compound %s for method %s" % ( compound, method ))
                    qval = None
                    zval = None
                    if method == PCGS or True:
                        ptype = particle['type']
                    else:
                        bases = particle['type'].split('_')
                        ptype = bases[0]
                        if len(bases) > 1:
                            if not bases[0].startswith('v'):
                                ptype = "v1"+bases[0]
                    try:
                        if method == "MBIS-S":
                            if compound in zeta[method] and ipart < len(zeta[method][compound]):
                                qval = float(data[method][compound][ipart]["value"])
                                zval = zeta[method][compound][ipart]["zeta"]
                        else:
                            for j in range(min(npart, len(data[method][compound]))):
                                if (ptype == data[method][compound][j]["type"] and
                                    particle["aindex"] == data[method][compound][j]["aindex"]):
                                    qval = float(data[method][compound][j]["value"])
                        if qval:
                            total[method] += qval
                            mystr.append( " & %s" % str(round(qval, 4)) )
                        else:
                            mystr.append( " & " )
                        if method in zeta and not method == "MBIS-S":
                            if compound in zeta[method] and ptype in zeta[method][compound]:
                                zval = zeta[method][compound][ptype]
                            else:
                                sys.exit("Cannot find zeta for particle %s in compound %s method %s" %
                                         ( particle, compound, method ) )

                        if zval:
                            mystr.append(" & %s" % ( str(round(zval, 2)) ))
                        elif method in zeta:
                            mystr.append(" & ")

                    except ValueError:
                        sys.exit("Don't know how to treat charge %s for type %s" % ( mypart["value"], mypart["type"] ))
                mystrings.append(mystr)
            
            # Now we have collected all atoms for all methods in a matrix
            # First, swap atom and shell for MBIS-S
            for j in range(0, len(mystrings)-1, 2):
                qatom  = mystrings[j][4]
                zatom  = mystrings[j][5]
                mystrings[j][4] = mystrings[j+1][4]
                mystrings[j][5] = mystrings[j+1][5]
                mystrings[j+1][4] = qatom
                mystrings[j+1][5] = zatom
            # Now print the lines
            for mystr in mystrings:
                for j in range(len(mystr)):
                    file.write(" %s" % mystr[j])
                file.write("\\\\\n")
                
            file.write("\\hline\n")
            file.write("Total")
            for method in data.keys():
                file.write(" & %g " % ( round(total[method],2) ) )
                if log_files[method]["ncol"] == 2:
                    file.write(" & ")
            file.write("\\\\\n")
            file.write("\\hline\n")
            file.write("\\end{tabular}\n")
            file.write("\\end{sidewaystable}\n")

    print("LaTeX file %s generated successfully!" % combined_output_file)


generate_logs()
data = extract_data_from_log()
save_data_as_latex(data)

