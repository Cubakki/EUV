from pymatgen.io.xyz import XYZ
from pymatgen.io.cif import CifParser
from pymatgen.core import Structure,Molecule
from bond_generate import bond_generate

molecule=Molecule.from_file("./data/881852.xyz")
print(bond_generate(molecule))