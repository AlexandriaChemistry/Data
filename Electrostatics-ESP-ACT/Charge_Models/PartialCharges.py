#!/usr/bin/env python

import glob, os, sys, json, xmltodict

debug = False

compounds_of_interest = [
    "ammonium",
    "methylammonium",
    "ethylammonium",
    "formate",
    "acetate",
    "guanidinium",
    "imidazolium",
    "water"
    ]

log_files = { "ESP":        { "label": "ESP", "ncol": 1 },
              "PC-elec":    { "label": "PC$_{e}$", "ncol": 1 },
              "PC-allelec": { "label": "PC$_{ei}$", "ncol": 1 },
              "MBIS-S":     { "label": "MBIS-S", "ncol": 2 },
              "PC+SV-elec": { "label": "PC+SV$_{e}$", "ncol": 2 },
              "PC+GV-elec": { "label": "PC+GV$_{e}$", "ncol": 2 },
              "PC+GS-elec": { "label": "PC+GS$_{e,i}$", "ncol": 2 },
              }

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
               "ow_s": [ "v1ow" ]
              }
    for key in newrep:
        if atype in newrep[key]:
            return key
    return atype

def get_zeta_mbiss(compound:str):
    zeta = []
    fn = f"../AntechamberGaussian/MBIS/{compound}/CCSD.json"
    if os.path.exists(fn):
        with open(fn, "r") as inf:
            data = json.load(inf)
        atoms = "atom_names"
        zzz   = "valence_widths"
        if atoms in data and zzz in data:
            for a in range(len(data[atoms])):
                zeta.append( { "atom": data[atoms][a], "zeta": 2*data[zzz][a] })
                zeta.append( { "atom": data[atoms][a]+"_s", "zeta": 0 })

def get_zeta(model:str, compound:str):
    zeta = []
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
                                                    idl   = entry[ident].lower()[:-2]
                                                    zeta.append({ "atom": idl, "zeta": float(value) })
                                                elif mytype == "BONDCORRECTIONS":
                                                    if not entry[ident] in bcc:
                                                        bcc[entry[ident]] = {}
                                                    bcc[entry[ident]][ptype] = value
                                                else:
                                                    if not entry[ident] in eem:
                                                        eem[entry[ident]] = {}
                                                    eem[entry[ident]][ptype] = value
                        
    return zeta, eem, bcc
        
def extract_data_from_log():
    data = {}
    gw   = "guanidinium#water"
    for file_type in log_files:
        log_file        = f"{file_type}_MP2.log"
        data[file_type] = {compound: [] for compound in compounds_of_interest}
        if not os.path.isfile(log_file):
            print(f"File does not exist: {log_file}")
            continue

        try:
            with open(log_file, 'r') as file:
                read_data = False
                current_compound = None
                core_shell_dict = {}

                guanidinium_hack = False
                for line in file:
                    line = line.strip()

                    if "Name:" in line:
                        for compound in compounds_of_interest:
                            c2 = ( "%s#%s" % ( compound, compound ) )
                            # Hack since we do not have homodimer data for guanidinium
                            if compound == "guanidinium" and f"Name: {gw}" in line:
                                guanidinium_hack = True
                                c2 = gw
                            if f"Name: {c2}" in line:
                                current_compound = compound
                                read_data = True
                                print("Will read charges for %s from %s" % ( compound, log_file ))
                                break
                        if current_compound:
                            continue  

                    if "EPOT" in line:
                        read_data = False
                        guanidinium_hack = False
                    elif read_data:
                        columns = line.split()
                        if len(columns) > 4 and columns[0].isdigit():
                            aindex     = 1
                            atom_type  = simplify(columns[2]) + ( "-%d" % aindex )
                            if not (guanidinium_hack and "w" in atom_type ):
                                acm_value  = columns[3]
                                for mypart in data[file_type][current_compound]:
                                    if mypart["type"] == atom_type:
                                        aindex   += 1
                                        atom_type = simplify(columns[2]) + ( "-%d" % aindex )
                                data[file_type][current_compound].append( { "type": atom_type, "value": acm_value } )
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
            file.write(r"\caption{Partial charges q (e) and screening widths $\\zeta$ (1/nm) for " + compound + " from ESP, MBIS-S and ACT models. First line, atom, second line shell or virtual site.}")
            file.write("\n")
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
            npart = len(data["PC+GS-elec"][compound])
            total = {}
            for method in data.keys():
                total[method] = 0
            # We use homodimers for reference and therefore only one half of the atoms is sufficient.
            if compound == "guanidinium":
                maxpart = npart
            else:
                maxpart = int(npart/2)
            for ipart in range(maxpart):
                particle = data["PC+GS-elec"][compound][ipart]
                mystr    = particle["type"].replace("_", "\\_")
                for method in data.keys():
                    if not compound in data[method]:
                        sys.exit("Could not find compound %s for method %s" % ( compound, method ))
                    try:
                        if method in zeta:
                            natom = len(data[method][compound])
                            fval = None
                            if ipart < natom:
                                fval = float(data[method][compound][ipart]["value"])
                            else:
                                sys.exit("Not enough atoms (%d, expected %d) in %s %s" % ( natom, maxpart, method, compound ) )
                            total[method] += fval
                            mystr += ( " & %s" % str(round(fval, 4)) )
                            zval = None
                            # Hardcoding stuff, sigh.
                            if method == "MBIS-S":
                                if ipart < len(zeta[method][compound]):
                                    zval = zeta[method][compound][ipart]["zeta"]
                            else:
                                zval = None
                                if method == "PC+GS-elec":
                                    ppp = particle['type'].split('-')[0]
                                else:
                                    ppp = "v1"+particle['type'].split('_')[0]
                                for j in range(len(zeta[method][compound])):
                                    if zeta[method][compound][j]['atom'] == ppp:
                                        zval = zeta[method][compound][j]['zeta']
                                        break

#                            elif particle in zeta[method][compound]:
#                                zval = zeta[method][compound][particle]
                            if zval:
                                mystr += (" & %s" % ( str(round(zval, 2)) ))
                            else:
                                mystr += " & - "
                        elif ipart % 2 == 0:
                            fval = float(data[method][compound][int(ipart/2)]["value"])
                            total[method] += fval
                            mystr += ( " & %s " %  str(round(fval, 4)) )
                        else:
                            mystr += " & "
                    except ValueError:
                        sys.exit("Don't know how to treat charge %s for type %s" % ( mypart["value"], mypart["type"] ))
                    
                file.write("%s\\\\\n" % mystr)
                
            file.write("\\hline\n")
            file.write("Total")
            for method in data.keys():
                file.write(" & \\multicolumn{%d}{c}{%g}" % ( log_files[method]["ncol"], round(total[method],2) ) )
            file.write("\\\\\n")
            file.write("\\hline\n")
            file.write("\\end{tabular}\n")
            file.write("\\end{sidewaystable}\n")

    print("LaTeX file generated successfully!")


data = extract_data_from_log()
save_data_as_latex(data)
