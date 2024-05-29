from pyxtal.symmetry import Group
from molecule import MolecularObject
import os
from readconfig import ReadConfig
from combination import WyckoffCombinations
from pyxtal.molecular_crystal import molecular_crystal
from pyxtal.symmetry import Wyckoff_position
from initialization import Initializer
from pyxtal import pyxtal
from ase.visualize import view
import random
import numpy as np
from pymatgen.symmetry.groups import SpaceGroup
from main import clean
from pyxtal.molecule import pyxtal_molecule
from routine import Routine
from chgnet.model.model import CHGNet
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from caltool import CalculateEnergy
import subprocess
from structure import structure_type_converter
from structure import distance_check
from combination import WyckoffCombinations
from analysis import vasp_ml_compare
from analysis import check_trajectory
from structure import Perturbation

config = ReadConfig()

sg = config.space_group()
mol_list = config.mol_list()
mol_number = config.mol_number()

x = WyckoffCombinations(mol_list)

for i in sg:
    yy = x.mol_ratio_comb_list_sg([1, 3], i, 3)

# per_dict = {
#     "d_lat": 0.0,
#     "d_coor": 0.2,
#     "d_rot": 1.0
# }
# poscar = Poscar.from_file(os.getcwd()+'/POSCAR_current')
# struc = poscar.structure

# perturb = Perturbation(struc)
# x = perturb.uniform(0.2, 2, 0.0)
# print(x)
# check_trajectory(struc, per_dict)


# pyxtal_sg = Group(148)
# comb = pyxtal_sg.list_wyckoff_combinations([3, 18])[0]
# print(comb)

# file = '/Users/qizhang/Desktop/paper_new/VASP_ML_compare/GeSeNa_compare.csv'
# title = "aaa"

# plt, df, count_change_sign_1, average = vasp_ml_compare(file, title, data_number=100)
# print(average)

# x = pyxtal(molecular=True)
# y = pyxtal(molecular=True)
# x.from_random(3, 89, mol_list, [3, 9], sites=[['1d', '1c', '1b'], ['8p', '1a']], factor=200.0)
# x = molecular_crystal(3, sg, mol_list, [3, 9], sites=[['2i', '1a'], ['2i', '1h', '1g', '1f', '1e', '1d', '1c', '1b']], factor=200.0)
# x = molecular_crystal(3, 148, mol_list, [3, 18], sites=[['3b'], ['18f']], factor=200.0)

# mol_crystal = molecular_crystal(dim=3, group=89, molecules=mol_list, numMols=[3, 9], factor=200.0, sites=[['2e', '1a'], ['2g', '2g', '2f', '1d', '1c', '1b']])
# print(mol_crystal.valid)

# print(mol_crystal)
# mol_crystal = molecular_crystal(dim=3, group=31, molecules=mol_list, numMols=[2, 12], factor=200.0, sites=[['2a'], ['4b', '4b', '4b']])
# x.from_random(3, 4, mol_list, [2, 12], sites=[['2a'], ['2a', '2a', '2a', '2a', '2a', '2a']], factor=1.0)
# x.from_random(3, 148, mol_list, [3, 18], sites=[['3b'], ['18f']], factor=5000.0)
# print(mol_crystal.valid)

# py = structure_type_converter(x, "pymatgen")

# poscar = Poscar.from_file("/Users/qizhang/Desktop/paper_new/PSLi/BESTgatheredPOSCARS_order")
# structure = poscar.structure
# cif = CifWriter(py, 0.1)
# cif.write_file("test.cif")

# y_py = structure_type_converter(y, "pymatgen")

# poscar = Poscar(y_py)
# poscar.write_file("CONTCAR_test")
# print(y)




# print(mol_list)
# x.from_random(3, sg, [mol_list[0]], [2], sites=[['2a']])
# x_pymatgen = x.to_pymatgen()
# y.from_seed(x_pymatgen, molecules=mol_list)
# y_pymatgen = y.to_pymatgen()
# z = pyxtal_molecule(mol_list[0], fix=True)
# print(y)
# print(x_pymatgen)
# print('-------------------------')
# print(y_pymatgen)
# structure = Structure.from_file(os.getcwd()+'/(d).cif')

# struc = CalculateEnergy(structure)

# for i in range(200):
#     chgnet = CHGNet.load()
#     # structure = Structure.from_file(os.getcwd()+'/POSCAR_sym148')
#     structure = Structure.from_file(os.getcwd()+'/(d).cif')
#     structure2 = Structure.from_file(os.getcwd()+'/(a).cif')
#     struc = CalculateEnergy(structure)
#     struc2 = CalculateEnergy(structure2)
#     energy = struc.mlp_CHGNet()[1]
#     energy2 = struc2.mlp_CHGNet()[1]
#     print(energy)
#     print(energy2)

# x = WyckoffCombinations(mol_list)
# y = x.mol_ratio_comb_list_sg(mol_number, 14, 3)

# x_ase.from_random(3, sg, mol_list, [3, 18])
# z_ase = x_ase.to_ase()
# view(z_ase)
# print(y)
# file = "/Users/qizhang/Desktop/phonon_check/POSCAR"
# poscar = Poscar.from_file(file)
# structure = poscar.structure
# analyzer = SpacegroupAnalyzer(structure, symprec=0.1)
# space_group = analyzer.get_space_group_symbol()
