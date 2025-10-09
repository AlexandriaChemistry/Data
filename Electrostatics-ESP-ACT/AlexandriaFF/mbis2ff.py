#!/usr/bin/env python3
import os, json, sys
import xml.etree.ElementTree as ET

BOHR = 0.0529177
json_dir = "../AntechamberGaussian/MBIS/"
xml_file = "P+S.xml"
output_xml = "P+S_updated.xml"

def set_val(var, val:float):
    var.set("value", str(val))
    var.set("minimum", str(val))
    var.set("maximum", str(val))
    var.set("mutability", "Fixed")

def average_data(atomtypes:list, data:list)->list:
    adict = {}
    for c in range(len(atomtypes)):
        key = atomtypes[c]
        if not key in adict:
            adict[key] = { "sum": data[c][0], "n": 1 }
        else:
            adict[key]["sum"] += data[c][0]
            adict[key]["n"] += 1
    alist = []
    for c in range(len(atomtypes)):
        key = atomtypes[c]
        val = adict[key]["sum"]/adict[key]["n"]
        alist.append([val])
    return alist

obtype_updates = {
    "guanidinium": [
        ("c2", "c2g"),
        ("n2", "n2g"),
        ("n2", "n2g"),
        ("n2", "n2g"),
        ("hn", "hng"),
        ("hn", "hng"),
        ("hn", "hng"),
        ("hn", "hng"),
        ("hn", "hng"),
        ("hn", "hng")
    ],
    "imidazolium": [
        ("c2", "c2i1"),
        ("c2", "c2i1"),
        ("n2", "n2i"),
        ("c2", "c2i2"),
        ("n2", "n2i"),
        ("h4", "h4i"),
        ("h4", "h4i"),
        ("hn", "hni"),
        ("h5", "h5i"),
        ("hn", "hni")
    ],
    "acetate": [
        ("c3", "c3a"),
        ("c2", "c2a"),
        ("o2", "o2a1"),
        ("o2", "o2a2"),
        ("hc", "hca1"),
        ("hc", "hca1"),
        ("hc", "hca1")
#        ("hc", "hca2"),
#        ("hc", "hca3")
    ],
    "ammonium": [
        ("hn", "hna1"),
        ("n4", "n4a"),
        ("hn", "hna1"),
        ("hn", "hna1"),
        ("hn", "hna1")
#        ("hn", "hna2"),
#        ("hn", "hna2")
    ],
    "ethylammonium": [
        ("n4", "n4e"),
        ("c3", "c3e1"),
        ("hn", "hne1"),
        ("hn", "hne1"),
        ("hn", "hne1"),
#        ("hn", "hne2"),
#        ("hn", "hne3"),
        ("h1", "h1e1"),
        ("h1", "h1e1"),
#        ("h1", "h1e2"),
        ("c3", "c3e2"),
        ("hc", "hce1"),
        ("hc", "hce1"),
        ("hc", "hce1"),
#        ("hc", "hce2"),
#        ("hc", "hce3"),
    ],
    "formate": [
        ("o2", "o2f"),
        ("o2", "o2f"),
        ("c2", "c2f"),
        ("h2", "h2f"),
    ],
    "methylammonium": [
        ("n4", "n4m"),
        ("c3", "c3m"),
        ("hn", "hnm1"),
        ("hn", "hnm1"),
        ("hn", "hnm1"),   
#        ("hn", "hnm2"),
#        ("hn", "hnm3"),   
        ("h1", "h1m1"),
        ("h1", "h1m1"),
        ("h1", "h1m1"),                     
#        ("h1", "h1m2"),
#        ("h1", "h1m3"),                     
    ],
    "propanoate": [
        ("o2", "o2p"),
        ("o2", "o2p"),
        ("c3", "c3p1"),
        ("c3", "c3p2"),
        ("c2", "c2p"),
        ("hc", "hcp1"),
        ("hc", "hcp1"),
        ("hc", "hcp2"),
        ("hc", "hcp2"),
        ("hc", "hcp2") 
    ],
    "butanoate": [
        ("c3", "c3b1"),
        ("o2", "o2b"),
        ("o2", "o2b"),
        ("c2", "c2b"),
        ("c3", "c3b2"),
        ("hc", "hcb1"),
        ("hc", "hcb1"),
        ("hc", "hcb2"),
        ("hc", "hcb2"),
        ("c3", "c3b3"),
        ("hc", "hcb3"),
        ("hc", "hcb3"),
        ("hc", "hcb3") 
    ],
    "water": [
        ("O", "ow"),
        ("H", "hw"),
        ("H", "hw"),
    ],
    "bromide": [("Br", "Br-")],
    "chloride": [("Cl", "Cl-")],
    "fluoride": [("F", "F-")],
    "potassium-ion": [("K", "K+")],
    "lithium-ion": [("Li", "Li+")],
    "sodium-ion": [("Na", "Na+")],
}

tree = ET.parse(xml_file)
root = tree.getroot()

coulomb_block = root.find(".//interaction[@type='COULOMB']")
if coulomb_block is None:
    raise RuntimeError("No COULOMB interaction block found in XML!")

particle_block = root.find(".//particletypes")
if particle_block is None:
    raise RuntimeError("No particletypes found in XML!")

for molname in obtype_updates.keys():
    json_path = os.path.join(json_dir, f"{molname}/CCSD.json")
    if not os.path.exists(json_path):
        print(f"no JSON for {molname}, skipping...")
        continue

    with open(json_path) as f:
        jdata = json.load(f)

    atomtypes = [new for _, new in obtype_updates[molname]]

    # Fetch widths
    val_widths = jdata.get("valence_widths", [])
    if not val_widths:
        print(f"no val_widths for {molname}, skipping...")
        continue

    # Now fetch charges
    charges = jdata.get("charges", [])
    val_charges = jdata.get("valence_charges", [])

    charges     = average_data(atomtypes, charges)
    val_charges = average_data(atomtypes, val_charges)
    val_widths  = average_data(atomtypes, val_widths)

    zeta_values = []
    for v in val_widths:
        zeta_values.append(1/(2*BOHR*float(v[0])))

    if len(atomtypes) != len(zeta_values):
        print(f"mismatch in number of atomtypes ({len(atomtypes)}) and zeta values ({len(zeta_values)}) for {molname}")
        continue

    print(f"{molname}: using {len(zeta_values)} zeta values")

    for atype, zval in zip(atomtypes, zeta_values):
        paramlist = coulomb_block.find(f".//parameterlist[@identifier='{atype}_z']")
        if paramlist is None:
            print(f"no <parameterlist> found for {atype}_z")
            continue

        param_elem = paramlist.find("parameter[@type='zeta']")
        if param_elem is None:
            print(f"no <parameter> element of type='zeta' for {atype}_z")
            continue

        zval_rounded = round(zval, 6)
        set_val(param_elem, zval_rounded)

    for c in range(len(charges)):
        qcore = charges[c][0] - val_charges[c][0]
        atype = atomtypes[c]
        # Change vtype
        vatype = "v1" + atype
        paramlist = particle_block.find(f".//particletype[@identifier='{vatype}']")
        if paramlist is None:
            sys.exit("Cannot find particletype for %s" % vatype)
        vcharge = paramlist.find(f".//parameter[@type='charge']")
        set_val(vcharge, qcore)
        # Change atom
        paramlist = particle_block.find(f".//particletype[@identifier='{atype}']")
        if paramlist is None:
            sys.exit("Cannot find particletype for %s" % atype)
        acharge = paramlist.find(f".//parameter[@type='charge']")
        set_val(acharge, val_charges[c][0])

tree.write(output_xml, xml_declaration=True, short_empty_elements=True)
print(f"updated XML written to {output_xml}")
