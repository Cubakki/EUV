from pymatgen.io.xyz import XYZ
from pymatgen.io.cif import CifParser
from pymatgen.core import Structure,Molecule
from bond_generate import bond_generate
from seperator import seperator

Seperator=seperator("./data/CONTCAR.xyz")
Seperator.run("Sn","C")