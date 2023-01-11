from pymatgen.io.cif import CifFile,CifParser
from pymatgen.io.xyz import XYZ

Parser=CifParser("./data/881852.cif")
structure=Parser.get_structures()
xyz=XYZ(structure)
xyz.write_file("./data/881852.xyz")


