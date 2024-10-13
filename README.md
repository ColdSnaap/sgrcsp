[![DOI](https://zenodo.org/badge/805911060.svg)](https://zenodo.org/doi/10.5281/zenodo.11838813)

[Link to the publication](https://doi.org/10.48550/arXiv.2407.21337) 

Sent to Journal of Physics: Condensed Matter and has passed peer review. The finalized version has been submitted.

<h1 align="center">Space Group Restricted Crystal Structure Prediction</h1>

</h4>
Space Group Restricted Crystal Structure Prediction (SGRCSP) is a specialized tool designed for predicting crystal structures based on a given chemical composition. This tool is particularly suited for structures featuring rigid bodies or molecules. One of the major differences between SGRCSP and other popular Crystal Structure Prediction (CSP) packages, such as USPEX or CALYPSO, is the ability of SGRCSP to restrict atoms (or rigid bodies) to designated Wyckoff positions. This restriction ensures that symmetry is preserved throughout the simulation.

## Prerequisites

### Python Version
Python 3.9+

### Required Packages
- `pymatgen`
- `ase`
- `pyxtal`
- `chgnet`

The `pyxtal` package will automatically install ASE version 3.18.0 and Pymatgen 2024.3.1 when you run:
```sh
pip install pyxtal
```
However, since the ASE package on PyPI has not been updated for two years, you need to install the latest ASE version from their GitLab repository to ensure compatibility with chgnet. You can do this by running:
```sh
pip install git+https://gitlab.com/ase/ase
```

### PyXtal modification
SGRCSP uses PyXtal to generate symmetry-restricted structures. However, some modifications and bug fixes are required for the PyXtal package to support our desired functionality. Please copy the Python scripts from the pyxtal_modify folder to the corresponding PyXtal package directory.

### Bond configuration
To generate molecular crystal structures, you need to set the bond lengths between atoms within the molecules. This configuration should be defined in the bonds.json file located in the database folder of the PyXtal package. Only approximate bond lengths are necessary.


## Usage
[Documentation](https://sgrcsp.readthedocs.io/en/latest/) is under development.

### Input files
- `input.txt`
- Molecule files: `MOL_1.xyz`, `MOL_2.xyz`
- First-principles calculation inputs

Please refer to the format of input files in the `examples/Input` folder. Currently, only VASP is supported as the first-principles calculation package.

### Combinations of Molecules/Atoms in a Space Group
Given the chemical composition of a structure in a space group, you can see how many possible compositions are available.

**Example `input.txt`:**
```sh
# Molecules in the structures are MOL_1.xyz and MOL_2.xyz
MolFiles = MOL_1.xyz MOL_2.xyz
# Ratio between these two molecules is 1:3
MolNumber = 1 3
# Maximum number of molecules can be 3:9
RatioMaxMultiplier = 3
# Space groups to search for possible combinations
SpaceGroup = 1 4 7 28-31
```
**Example Python Script:**
```python
from readconfig import ReadConfig
from combination import WyckoffCombinations

config = ReadConfig()
sg = config.space_group()
mol_list = config.mol_list()
mol_number = config.mol_number()
z = config.max_multiplier()

comb = WyckoffCombinations(mol_list)
for i in sg:
    comb_number = comb.mol_ratio_comb_list_sg(mol_number, i, z)
```
**Example Output:**
```python
P1 (1), Z = 1, checking: 1/1, valid: 1, sites: [['1a'], ['1a', '1a', '1a']]                      
P1 (1), Z = 2, checking: 1/1, valid: 1, sites: [['1a', '1a'], ['1a', '1a', '1a', '1a', '1a', '1a']]                      
P1 (1), Z = 3, checking: 1/1, valid: 1, sites: [['1a', '1a', '1a'], ['1a', '1a', '1a', '1a', '1a', '1a', '1a', '1a', '1a']]                      
P21 (4), Z = 2, checking: 1/1, valid: 1, sites: [['2a'], ['2a', '2a', '2a']]                      
Pc (7), Z = 2, checking: 1/1, valid: 1, sites: [['2a'], ['2a', '2a', '2a']]                      
Pma2 (28), Z = 2, checking: 39/39, valid: 39, sites: [['2c'], ['4d', '2c']]                            
Pnc2 (30), Z = 2, checking: 12/12, valid: 12, sites: [['2b'], ['4c', '2b']]                            
Pmn21 (31), Z = 2, checking: 2/2, valid: 2, sites: [['2a'], ['4b', '2a']]  
```

### Space group restricted crystal structure prediction
1.	Run __init__.py to generate a simulation folder for each valid combination.
2.	Run main.py in each folder to start the simulation.
