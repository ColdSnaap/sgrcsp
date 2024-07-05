import os
from collections import Counter
from readconfig import ReadConfig
from version import __version__
from initialization import Initializer
from combination import WyckoffCombinations

config = ReadConfig()
sg_lsit = config.space_group()
sg_mol = config.mol_list()

print("Version:", __version__)
print("Spacegroup:", config.config["SpaceGroup"])

for i, mol in enumerate(sg_mol):
    elements = mol.labels
    element_counts = Counter(elements)
    formula = ''.join(f"{element}{count if count > 1 else ''}" for element, count in element_counts.items())
    print(f"Molecule{i+1}: {formula}")

print("\n")

print("Searching...")
# remove mob.txt if exist
comb_result = os.getcwd()+"/comb.txt"
if os.path.exists(comb_result):
    os.remove(comb_result)

for i, sg in enumerate(sg_lsit):
    mkdir = Initializer(sg).make_comb_dir()