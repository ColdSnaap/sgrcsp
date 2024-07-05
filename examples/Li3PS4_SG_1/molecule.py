import numpy as np
from pymatgen.core.structure import IMolecule

class MolecularObject:
    """
    Generate various types of molecular objects from a file.
    pymatgen: pymatgen molecular object 
    """
    def __init__(self, mol_file: str) -> None:
        self.mol_file = mol_file
        self.file_type = self.mol_file.split('.')[-1]

    def read_molecular_file(self):
        """
        xyz: file with xyz cartesian coordinates of atoms
            - Add fractional or cartesian coordiante check
        """
        if self.file_type == 'xyz':
            try:
                with open(self.mol_file, 'r') as file:
                    species = file.readline().strip().split()
                    if len(species) == 1:
                        coords_array = np.loadtxt(self.mol_file, skiprows=1).reshape(1,3)
                    else:
                        coords_array = np.loadtxt(self.mol_file, skiprows=1)
                return species, coords_array
            except Exception as e:
                raise Exception(f"Failed to read molecular file: {e}")

    def pymatgen(self):
        species, coords = self.read_molecular_file()
        if species and coords is not None:
            molecular = IMolecule(species, coords)
            return molecular
        else:
            raise ValueError("No data found in file to create molecular object.")

