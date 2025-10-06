#!/usr/bin/env python

import xml.etree.cElementTree as ET

tree = ET.parse("P+S.xml")
ac   = "alexandria_chemistry"
# Gets root tag on XML
root = tree.getroot()
for ind,val in enumerate(root):
    # Checking the tags for each CD element
    for elem in val.iter():
        # If tag is year and value of tag is 1990
        if elem.tag == 'interaction' and elem.get('function') == "COULOMB_SLATER":
            for plist in elem.findall('parameterlist'):
                for param in plist.findall('parameter'):
                    for change in [ 'value', 'maximum' ]:
                        param.set(change, str(2*float(param.get(change))))
tree.write("P+S-hacked.xml", xml_declaration=True, encoding='utf-8')
