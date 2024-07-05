import os
from readconfig import ReadConfig
from pymatgen.io.vasp import Poscar
from structure import apply_perturbation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter

file = "/Users/qizhang/Desktop/phonon_check/structure_check/POSCAR3"
poscar = Poscar.from_file(file)
structure = poscar.structure

sga = SpacegroupAnalyzer(structure, symprec=0.1)
print(sga.get_space_group_number())

refined_structure = sga.get_refined_structure()


cif_writer = CifWriter(refined_structure, symprec=0.1)
cif_writer.write_file("/Users/qizhang/Desktop/phonon_check/structure_check//poscar3.cif")