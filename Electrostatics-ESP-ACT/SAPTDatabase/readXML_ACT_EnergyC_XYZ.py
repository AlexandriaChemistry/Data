#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import os
import csv

xmlfile = "../AlexandriaFF/sapt-esp.xml"

dimers = [
    "water#water", "water#bromide", "water#lithium", "water#potassium", "potassium#chloride",
    "acetate#lithium", "acetate#sodium", "acetate#potassium", "acetate#water", "acetate#fluoride",
    "acetate#chloride", "methylammonium#methylammonium", "methylammonium#fluoride", "methylammonium#chloride",
    "methylammonium#bromide", "methylammonium#water", "methylammonium#lithium", "methylammonium#sodium",
    "ammonium#acetate", "ammonium#ammonium", "ammonium#bromide", "ammonium#water", "ammonium#lithium",
    "ammonium#sodium", "ethylammonium#bromide", "ethylammonium#chloride", "ethylammonium#ethylammonium",
    "ethylammonium#water", "formate#formate", "formate#potassium", "formate#water", "sodium#sodium",
    "potassium#potassium", "chloride#chloride", "bromide#bromide", "sodium#chloride", "chloride#fluoride",
    "formate#chloride", "potassium#bromide", "ethylammonium#potassium", "sodium#fluoride", "lithium#bromide",
    "lithium#sodium", "bromide#fluoride", "bromide#chloride", "water#chloride", "water#fluoride",
    "water#sodium", "ethylammonium#methylammonium", "sodium#potassium", "formate#sodium", "formate#lithium",
    "acetate#acetate", "acetate#bromide", "methylammonium#formate", "methylammonium#acetate", "ammonium#chloride",
    "ammonium#fluoride", "ammonium#formate", "ethylammonium#fluoride", "ethylammonium#formate", "lithium#lithium",
    "fluoride#fluoride", "lithium#fluoride",  "ammonium#potassium", "ethylammonium#lithium",
    "ethylammonium#sodium", "formate#fluoride", "formate#bromide", "methylammonium#potassium", "potassium#fluoride",
    "sodium#bromide", "lithium#chloride", "lithium#potassium", "ammonium#ethylammonium", "formate#acetate",
    "ammonium#methylammonium"
]

def writeXYZ(atom_data, xyzfilename):
    with open(xyzfilename, 'w') as f:
        f.write(f"{len(atom_data)}\n")
        f.write("Coordinates from XML\n")
        for atom in atom_data:
            f.write(f"{atom['name']} {atom['x']} {atom['y']} {atom['z']}\n")

def parseXMLandgiveXYZ(xmlfile, xyz_dir):
    tree = ET.parse(xmlfile)
    root = tree.getroot()


    if not os.path.exists(xyz_dir):
        os.makedirs(xyz_dir)

    for dimer in dimers:
        molecule = root.find(f".//molecule[@molname='{dimer}']")

        if molecule is not None:
            for experiment in molecule.findall(".//experiment"):
                datafile = experiment.attrib["datafile"]


                dimer_dir = os.path.join(xyz_dir, dimer)
                if not os.path.exists(dimer_dir):
                    os.makedirs(dimer_dir)

                atom_data = []
                stop_collecting = False

                for elem in experiment.iter():
                    if stop_collecting:
                        break

                    if elem.tag == "atom":
                        atom_data.append({
                            'name': elem.get('name'),
                            'x': elem.find('x').text,
                            'y': elem.find('y').text,
                            'z': elem.find('z').text
                        })

                    if elem.tag == 'atom' and elem.tail == '</atom>\n  </experiment>':
                        stop_collecting = True


                if atom_data:
                    xyzfilename = os.path.join(dimer_dir, f"{datafile}.xyz")
                    writeXYZ(atom_data, xyzfilename)

def extract(xmlfile):
    tree = ET.parse(xmlfile)
    root = tree.getroot()

    with open('energies-SC-Water-Ions.csv', mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='|')

        energycomponents = [
            "Disp2 (CCD)", "Disp20", "Disp21", "Disp22 (S) (CCD)", "Disp22 (SDQ)", "Disp22 (T)", "Disp22 (T) (CCD)",
            "Dispersion", "Electrostatics", "Elst10,r", "Elst12,r", "Est. Disp22 (T)", "Exch-Disp20", "Exch-Ind20,r",
            "Exch-Ind22", "Exch10", "Exch10(S^2)", "Exch11(S^2)", "Exch12(S^2)", "Exchange", "Ind20,r", "Induction",
            "InteractionEnergy", "delta HF,r (2)", "delta MP2,r (2)"
        ]

        header = ["Molecule", "Datafile"] + energycomponents
        writer.writerow(header)

        for dimer in dimers:
            molecule = root.find(f".//molecule[@molname='{dimer}']")

            if molecule is not None:
                molname = molecule.attrib['molname']

                for experiment in molecule.findall(".//experiment"):
                    datafile = experiment.attrib["datafile"]

                    energy_values = []

                    for energy_type in energycomponents:
                        energy_tag = experiment.find(f".//energy[@type='{energy_type}']/average")
                        if energy_tag is not None:
                            energy_values.append(energy_tag.text)
                        else:
                            energy_values.append("N/A")

                    writer.writerow([molname, datafile] + energy_values)

    print("Yeeeeeeeees! data extraction complete! Results are saved to energies-SC-Water-Ions.csv.")


xyzdir = "XYZ"
parseXMLandgiveXYZ(xmlfile, xyzdir)
print("Hooray! Now you have all your desired XYZ files!")
extract(xmlfile)

