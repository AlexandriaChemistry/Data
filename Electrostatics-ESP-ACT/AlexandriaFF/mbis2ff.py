#!/usr/bin/env python3
import os
import json
import xml.etree.ElementTree as ET

json_dir = "."
xml_file = "P+S.xml"
output_xml = "P+S_updated.xml"

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
        ("hn", "hng"),
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
        ("hn", "hni"),        
    ],
    "acetate": [
        ("c3", "c3a"),
        ("c2", "c2a"),
        ("o2", "o2a1"),
        ("o2", "o2a2"),
        ("hc", "hca1"),
        ("hc", "hca2"),
        ("hc", "hca3"),
    ],
    "ammonium": [
        ("hn", "hna1"),
        ("n4", "n4a"),
        ("hn", "hna1"),
        ("hn", "hna2"),
        ("hn", "hna2"),
    ],
    "ethylammonium": [
        ("n4", "n4e"),
        ("c3", "c3e1"),
        ("hn", "hne1"),
        ("hn", "hne2"),
        ("hn", "hne3"),
        ("h1", "h1e1"),
        ("h1", "h1e2"),
        ("c3", "c3e2"),
        ("hc", "hce1"),
        ("hc", "hce2"),
        ("hc", "hce3"),
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
        ("hn", "hnm2"),
        ("hn", "hnm3"),   
        ("h1", "h1m1"),
        ("h1", "h1m2"),
        ("h1", "h1m3"),                     
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

for molname in obtype_updates.keys():
    json_path = os.path.join(json_dir, f"{molname}_mbis_ps.json")
    if not os.path.exists(json_path):
        print(f"no JSON for {molname}, skipping...")
        continue

    with open(json_path) as f:
        jdata = json.load(f)

    sigma_inv_nm = jdata.get("sigma_inv_nm", [])
    if not sigma_inv_nm:
        print(f"no sigma_inv_nm for {molname}, skipping...")
        continue

    zeta_values = [v/2 for v in sigma_inv_nm]

    atomtypes = [new for _, new in obtype_updates[molname]]

    if len(atomtypes) != len(zeta_values):
        print(f"mismatch in number of atomtypes ({len(atomtypes)}) and zeta values ({len(zeta_values)}) for {molname}")
        continue

    print(f"{molname}: using {len(zeta_values)} zeta values (divided by 2)")

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
        param_elem.set("value", str(zval_rounded))
        param_elem.set("minimum", str(zval_rounded))
        param_elem.set("maximum", str(zval_rounded))

tree.write(output_xml, xml_declaration=True, short_empty_elements=True)
print(f"updated XML written to {output_xml}")
