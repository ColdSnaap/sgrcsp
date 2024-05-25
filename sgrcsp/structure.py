import os
import copy
import sys
from itertools import combinations
from pathlib import Path
from readconfig import ReadConfig
from pyxtal import pyxtal
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifParser


def structure_type_converter(structure, target_type, molecule=True):
    """
    Strcuture type converter.
    If structure is a path such as POSCAR or cif file,
    this function will read and store the structure as a pymatgen structure object first.
    Support types:
        Pymatgen <- File path, pymatgen, pyxtal
        Pyxtal <- File path, Pyxtal, pymatgen

    Arg:
        structure: various structure types (check support types in the description)
        target_type: targeted structure type
        molecule: convert to molecule pyxtal structure
            (molecule should always be True since single atoms are treated
            as molecules with one atom in this package)
    """
    
    if isinstance(structure, str):
        # Create a Path object from the file path
        path = Path(structure)
        if path.exists():
            # Check if the file has an extension
            if path.suffix:
                # Return the extension without the dot
                file_type = path.suffix[1:]
            else:
                # Return the file name if there is no extension
                file_type = path.name
        else:
            raise NameError("File not exist")
        
        # Read structure file and store it as pymatgen structure
        if file_type == "cif":
            parser = CifParser(structure)
            structure = parser.get_structures()[0]
        elif file_type == "POSCAR":
            poscar = Poscar.from_file(structure)
            structure = poscar.structure
        
        if target_type == "pyxtal":
            pyxtal_struc = pyxtal(molecular=molecule)
            if molecule:
                config = ReadConfig()
                mol_list = config.mol_list()
                pyxtal_struc.from_seed(structure, molecules=mol_list)
                return pyxtal_struc
            else:
                raise TypeError(f"molecule={molecule} not support yet")
        elif target_type == "pymatgen":
            return structure
                
    if target_type == "pymatgen":
        if isinstance(structure, Structure):
            return structure
        elif isinstance(structure, pyxtal):
            structure = structure.to_pymatgen()
            return structure

    if target_type == "pyxtal":
        if isinstance(structure, Structure):
            pyxtal_struc = pyxtal(molecular=molecule)
            if molecule:
                config = ReadConfig()
                mol_list = config.mol_list()
                pyxtal_struc.from_seed(structure, molecules=mol_list)
                return pyxtal_struc
            else:
                raise TypeError(f"molecule={molecule} not support yet")
        elif isinstance(structure, pyxtal):
            return structure
    
    else:
        raise NameError(f"Does not support target_type={target_type} yet")


def distance_check(structure) -> float:
    """
    Return the minimum distance between atoms in a structure.
    """
    
    structure = structure_type_converter(structure, "pymatgen")

    min_distance = float('inf')  # Initialize with infinity

    # Generate all unique pairs of sites
    for site1, site2 in combinations(structure.sites, 2):
        # Calculate distance between the pair of sites directly
        distance = site1.distance(site2)
        # Update the minimum distance found
        if distance < min_distance:
            min_distance = distance

    return min_distance


def apply_perturbation(structure, d_coor, d_rot, d_lat, trial_max=1000):

    struc_pyxtal = structure_type_converter(structure, "pyxtal")
    struc_trial = copy.deepcopy(struc_pyxtal)
    config = ReadConfig()
    min_dis_allow = config.min_dis()

    min_dis_struc = 0.001

    print("Applying perturbation")
    print(f"d_lat, d_coor, d_rot: {d_lat}, {d_coor}, {d_rot}")
    sys.stdout.flush()
    count = 0
    
    while min_dis_struc < min_dis_allow:
        count += 1
        struc_trial.apply_perturbation(d_lat, d_coor, d_rot)
        
        struc_py = structure_type_converter(struc_trial, "pymatgen")
        min_dis_struc = distance_check(struc_py)

        if count == trial_max:
            print(f"Could not find good structure (dis_allowed={min_dis_allow}, dis_found={min_dis_struc}, trial_max={trial_max})")
            print(f"Appling smaller perturbations")
            sys.stdout.flush()
            d_lat *= 0.9
            d_coor *= 0.9
            d_rot *= 0.9
            print(f"d_lat, d_coor, d_rot: {d_lat}, {d_coor}, {d_rot}")

            count = 0
    
    print(f"Trail: {count}")
    sys.stdout.flush()

    return struc_trial


def structure_generation(
    sg,
    mol_list,
    mol_number,
    sites,
    factor=1.0
) -> pyxtal:
    print("Generating structure ...")
    sys.stdout.flush()
    while True:
        try:
            structure = pyxtal(molecular=True)
            structure.from_random(3, sg, mol_list, mol_number, sites=sites, factor=factor)
            break
        except RuntimeError as e:
            print(f"RuntimeError: {e}.")
        finally:
            factor += 0.05
            print(f"Increasing volume, volume factor: {factor}")
            sys.stdout.flush()
    
    print("Structure generated\n")
    return structure


# def relaxed_structure_with_intact_molecule(
#     structure,
#     caltool,
#     max_trial=10
# ) -> pyxtal:
#     """
#     Check if the molecules are intact after structure relaxation.
#     If not generate a new one structrue and relax it again until reach max_trial.
    
#     Args:
#         structure: structure object
#         (To see supported type please check structure_type_converter function)
#         caltool: structure relaxation tool
#         max_trial: Maximum attempts to generate structure with intact molecules
    
#     Return:
#         pyxtal structure object
#     """

#     config = ReadConfig()
#     sub_command = config.sub_command
#     cal = CalculateEnergy(structure)
#     relax_struc = cal.energy_calculate(caltool, True, sub_command)[0]
#     structure = structure_type_converter(relax_struc, "pyxtal")
#     pass