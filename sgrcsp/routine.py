import os
import shutil
import sys
import ast
from collections import Counter
from optimization import GlobalOptimize
from readconfig import ReadConfig
from caltool import CalculateEnergy
from pyxtal import pyxtal
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp.inputs import Poscar
from pyxtal.msg import ReadSeedError
from structure import structure_type_converter
from structure import structure_generation
from structure import apply_perturbation


def modify_update(file_path, key, new_value):
    # Define the key to look for in the format it appears in the file (e.g., "SpaceGroup =")
    key_with_equal = f"{key} ="

    # Read the existing content of the file
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print("File not found.")
        return

    # Modify the line containing the key
    with open(file_path, 'w') as file:
        for line in lines:
            # Check if the line contains the key and modify it
            if line.strip().startswith(key_with_equal):
                line = f"{key_with_equal} {new_value}\n"
            file.write(line)


def continue_log(
    step_current,
    accept_list,
    current_state,
    count,
    best_energy,
    best_state
):
    continue_path = os.getcwd() + "/Result/Continue"
    continue_file = continue_path + "/continue"
    if not os.path.isdir(continue_path):
        os.makedirs(continue_path)
    if not os.path.isfile(continue_file):
        with open(continue_file, 'w') as cf:
            cf.write("#Continue config")
        with open(continue_file, 'a') as file:
            file.write(f'\nstep_current = {step_current}')
            file.write(f'\naccept_list = {accept_list}')
            file.write(f'\ncount = {count}')
            file.write(f'\nbest_energy = {best_energy}')
    else:
        modify_update(continue_file, "step_current", step_current)
        modify_update(continue_file, "accept_list", accept_list)
        modify_update(continue_file, "count", count)
        modify_update(continue_file, "best_energy", best_energy)

    current_state = structure_type_converter(current_state, "pymatgen")
    best_state = structure_type_converter(best_state, "pymatgen")
    poscar_current = Poscar(current_state)
    poscar_current.write_file(continue_path+"/POSCAR_current")
    poscar_best = Poscar(best_state)
    poscar_best.write_file(continue_path+"/POSCAR_best")
    

def read_continue_config():
    continue_path = os.getcwd() + "/Result/Continue"
    continue_file = continue_path + "/continue"
    config = {}
    try:
        with open(continue_file, 'r') as file:
            for line in file:
                # Ignore empty lines and lines starting with '#'
                if not line or line.startswith('#'):
                    continue
                # Ignore empty lines and lines without '='
                if "=" not in line:
                    continue
                key, value = line.strip().split('=', 1)
                config[key.strip()] = value.strip()
    except FileNotFoundError:
        print("Configuration file not found.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
    
    config["step_current"] = int(config["step_current"])
    config["count"] = int(config["count"])
    config["best_energy"] = float(config["best_energy"])
    config["accept_list"] = ast.literal_eval(config["accept_list"])
    current_state = Poscar.from_file(continue_path+"/POSCAR_current")
    best_state = Poscar.from_file(continue_path+"/POSCAR_best")
    
    return config["step_current"], config["accept_list"], config["count"], config["best_energy"], \
        current_state, best_state


class Routine:
    """
    Determined routine for gloable optimization using the input.txt from each combination.
    SpaceGroup = int()
    """
    def __init__(self) -> None:
        self.config = ReadConfig()
        self.routine = self.config.routine_number()
        self.sites = self.config.mol_sites()
        self.sg = self.config.space_group()[0]
        self.mol_list = self.config.mol_list()
        self.mol_number = self.config.mol_number()
        self.max_step = self.config.max_step()
        self.inter_loop = self.config.internal_loop()
        self.job_contine = self.config.job_continue()
        self.sub_command = self.config.sub_command()
        self.factor = self.config.volume()
        self.caltool = self.config.cal_tool()
        self.initial_relax = self.config.initial_relax()

        self.loop_number = int(self.max_step/self.inter_loop)
        # input information
        print(f"Routine: {self.routine}")
        print(f"Spacegroup: {self.sg}")
        for i, mol in enumerate(self.mol_list):
            elements = mol.labels
            element_counts = Counter(elements)
            formula = ''.join(f"{element}{count if count > 1 else ''}" for element, count in element_counts.items())
            print(f"Molecule{i+1}: {formula}")
        print(f"Sites: {self.sites}")
        print(f"Number: {self.mol_number}")
        print(f"Volume: {self.factor}")
        print(f"LoopNumber: {self.inter_loop}")
        print(f"MaxStep: {self.max_step}\n")
        sys.stdout.flush()
        # generate initial structure
        if self.initial_relax:
            while True:
                try:
                    self.structure = structure_generation(
                        self.sg,
                        self.mol_list,
                        self.mol_number,
                        self.sites,
                        factor=self.factor
                    )
                    print("Relaxing the initial structure ...")
                    sys.stdout.flush()
                    cal = CalculateEnergy(self.structure)
                    cal_struc = cal.energy_calculate("VASP", True, self.sub_command)[0]
                    self.structure = structure_type_converter(cal_struc, "pyxtal")
                    print("Initial structure generated")
                    sys.stdout.flush()
                    break
                except ReadSeedError:
                    print("Molecule not intact, trying to generate a new structure to relax")
        else:
            self.structure = structure_generation(
                self.sg,
                self.mol_list,
                self.mol_number,
                self.sites,
                factor=self.factor
            )                 

    def log(
        self,
        step: int,
        structure,
        energy: float,
        trajectory: bool = False,
    ) -> None:
        if isinstance(structure, pyxtal):
            structure = structure.to_pymatgen()
        
        if not isinstance(structure, Structure):
            print(type(structure))
            raise TypeError("Structure is not pymatgen structure type")
        
        result_dir = os.getcwd() + "/Result"
        struc_log = os.getcwd() + "/Result/BestStrucsList"
        struc_file = os.getcwd() + "/Result/BestStrucs"
        
        sga = SpacegroupAnalyzer(structure, symprec=0.1)
        refined_structure = sga.get_refined_structure()
        space_group_number = sga.get_space_group_number()
        space_group_symbol = sga.get_space_group_symbol()
        
        if not os.path.isdir(result_dir):
            os.mkdir(result_dir)
            # raise NameError("Result directory not exist")
        if not os.path.isdir(struc_file):
            os.makedirs(struc_file)

        name_en = abs(energy)
        if not os.path.isfile(f"{struc_file}/{name_en}.cif"):
            # refined_structure.to(f'{struc_file}/{name_en}.cif', fmt="cif")
            cif_writer = CifWriter(refined_structure, symprec=0.1)
            cif_writer.write_file(f'{struc_file}/{name_en}.cif')
        if not os.path.isfile(struc_log):
            with open(struc_log, 'w') as log_write:
                log_write.write("Steps, Space Group, Energy/Atom\n")
                log_write.write(f'{step}, {space_group_symbol}({space_group_number}), {energy}\n')
        else:
            with open(struc_log, 'a') as log_write:
                log_write.write(f'{step}, {space_group_symbol}({space_group_number}), {energy}\n')


    # def routine1(self):
    #     """
    #     Routine1:
    #     1. Relax the initial structure with vasp.
    #     2. Random sampling for a certain step.
    #     3. Relax and store the structure after random sampling.
    #     4. Start from relaxed structure and give it a small perturbation
    #     5. Repeat 2-3 for a certain step.
    #     6. Order and store relaxed structures.
    #     """
    #     # Relax the initial structure with vasp.
    #     sub_command = self.config.sub_command()
    #     relax_intial = self.energy_calculate(self.structure, 'VASP', True, sub_command)
    #     if relax_intial[0] is False:
    #         raise TypeError("Structure are not relaxed")
    #     else:
    #         structure = relax_intial[0]

    #     # Random sampling for a certain step
    #     struc_pyxtal1 = pyxtal(molecular=True)
    #     struc_pyxtal = struc_pyxtal1.from_seed(structure, molecules=self.mol_list)
    #     optimizer = GlobalOptimize(struc_pyxtal)
    #     for i in range(int(self.max_step/self.inter_loop)):
    #         optimizer = GlobalOptimize(struc_pyxtal)
    #         if i == 0:
    #             scale_dict = optimizer.sampling_scale(struc_pyxtal)
    #             output = optimizer.simulated_annealing(scale_dict=scale_dict)
    #             print(f'scale_list = {scale_dict}')
    #             step_current, temp_current, accept_list, current_state = output
    #         else:
    #             scale_dict = optimizer.sampling_scale(current_state)
    #             output = optimizer.simulated_annealing(step_current, temp_current, accept_list, scale_dict)
    #             print(f'scale_list = {scale_dict}')
    #             step_current, temp_current, accept_list, current_state = output

    #         # Relax and store the structure after random sampling
    #         # Store file and information of relaxed to the result direcotry 
    #         struc_relaxed = self.energy_calculate(self.structure, 'VASP', True, sub_command)
    #         self.log(step_current, struc_relaxed, struc_relaxed[2])
    

    def routine2(self):
        """
        Routine2:
        Full simulated annealing for debug.
        """
        print("Simulated annealing without structure relaxation")
        print(self.job_contine)
        if not self.job_contine:
            # clean Log file
            result_dir = os.getcwd() + "/Result"
            if os.path.isfile(result_dir+"/Log"):
                os.remove(result_dir+"/Log")
            if os.path.isfile(result_dir+"/BestStrucsList"):
                os.remove(result_dir+"/BestStrucsList")
            if os.path.isdir(result_dir+"/BestStrucs"):
                shutil.rmtree(result_dir+"/BestStrucs")
            
            step_current = 0
            accept_list = None
            current_state = self.structure
            count = None
            best_energy = None
            best_state = None
        else:
            pass

        for step in range(self.loop_number):
            optimizer = GlobalOptimize(current_state)
            output = optimizer.simulated_annealing(step_current, accept_list, count, best_energy, best_state, big_loop=step, trajectory=True)
            step_current, temp_current, accept_list, current_state, count, best_energy, best_state = output

            # Store best structure
            self.log(step_current, best_state, best_energy)
            print(f"step:{step_current}, temp:{temp_current}")
    
   
    def routine3(self):
        """
        Routine3:
        relax using VASP, SA using ML.
        """        
        if not self.job_contine:
            # clean Log file
            result_dir = os.getcwd() + "/Result"
            if os.path.isfile(result_dir+"/Log"):
                os.remove(result_dir+"/Log")
            if os.path.isfile(result_dir+"/BestStrucsList"):
                os.remove(result_dir+"/BestStrucsList")
            if os.path.isdir(result_dir+"/BestStrucs"):
                shutil.rmtree(result_dir+"/BestStrucs")
            
            step_current = 0
            accept_list = None
            current_state = self.structure
            count = None
            best_energy = None
            best_state = None
            rate_scale = None
        else:
            step_current, accept_list, count, best_energy, current_state, best_state = \
                read_continue_config()
        
        per_dict = {
            "d_lat": 0.0,
            "d_coor": 0.1,
            "d_rot": 1.0
        }
        

        for step in range(self.loop_number):
            
            while True:
                try:
                    current_state = structure_type_converter(current_state, "pymatgen")
                    cal1 = CalculateEnergy(current_state)
                    cal = cal1.energy_calculate("VASP", True, self.sub_command)
                    current_state, relax_energy, relax_energy_atom = cal
                    current_state = structure_type_converter(current_state, "pyxtal")
                    break
                except ReadSeedError:
                    print("Molecule not intact, trying to generate a new structure to relax")
                    sys.stdout.flush()
                    current_state = apply_perturbation(current_state, 0.2, 1.0, 0.0)
            
            self.log(step_current, current_state, relax_energy)     

            # small perturbation
            current_state = apply_perturbation(current_state, 0.1, 1.0, 0.0)

            optimizer = GlobalOptimize(current_state)
            output = optimizer.simulated_annealing(
                step_current,
                accept_list,
                count,
                best_energy,
                best_state,
                big_loop=step,
                rate_scale=rate_scale,
                perturb_dict = per_dict
            )
            step_current, temp_current, accept_list, current_state, count, best_energy, best_state, rate_scale = output

            continue_log(step_current, accept_list, current_state, count, best_energy, best_state)

            if step == self.loop_number - 1:
                print(f"Reach max loop: {step}")
            else:
                print(f"step:{step_current}, temp:{temp_current}")
    

    def routine4(self):
        """
        Routine3:
        relax using VASP, SA using VASP.
        """        
        if not self.job_contine:
            # clean Log file
            result_dir = os.getcwd() + "/Result"
            if os.path.isfile(result_dir+"/Log"):
                os.remove(result_dir+"/Log")
            if os.path.isfile(result_dir+"/BestStrucsList"):
                os.remove(result_dir+"/BestStrucsList")
            if os.path.isdir(result_dir+"/BestStrucs"):
                shutil.rmtree(result_dir+"/BestStrucs")
            
            step_current = 0
            accept_list = None
            current_state = self.structure
            count = None
            best_energy = None
            best_state = None
        else:
            step_current, accept_list, count, best_energy, current_state, best_state = \
                read_continue_config()
        
        for step in range(self.loop_number):
            # structure relax
            current_state = structure_type_converter(current_state, "pymatgen")
            cal1 = CalculateEnergy(current_state)
            cal = cal1.energy_calculate("VASP", True, self.sub_command)
            current_state, relax_energy, relax_energy_atom = cal
            self.log(step_current, current_state, relax_energy)     

            optimizer = GlobalOptimize(current_state)
            output = optimizer.simulated_annealing(step_current, accept_list, count, best_energy, best_state, big_loop=step)
            step_current, temp_current, accept_list, current_state, count, best_energy, best_state = output

            continue_log(step_current, accept_list, current_state, count, best_energy, best_state)

            print(f"step:{step_current}, temp:{temp_current}")