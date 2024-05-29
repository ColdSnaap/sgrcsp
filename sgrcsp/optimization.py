import os
import math
import random
import copy
import sys
import shutil
import numpy as np
from caltool import CalculateEnergy
from readconfig import ReadConfig
from pyxtal import pyxtal
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Poscar
from structure import structure_type_converter
from structure import apply_perturbation

class GlobalOptimize:

    def __init__(self, structure) -> None:
        
        self.structure = structure_type_converter(structure, "pyxtal", True)
        if not isinstance(self.structure, pyxtal):
            raise TypeError("Structure should be pyxtal structure object")

        self.config = ReadConfig()
    

    def cooling_rate(self, temp_initial: float, function, step, temp_final = None):
        if function == 'constant':
            temp = temp_initial
        elif function == 'fast':
            temp = temp_initial * 0.99**(step)
        return round(temp, 3)


    def random_sampling(self, temp, function, constant = None, sample_number = 1):
        if function == 'constant':
            result = constant
        elif function == 'normal':
            if sample_number == 1:
                # result is a float
                result = abs(float(np.random.normal(0, temp, sample_number)[0]))
            else:
                # result is a list
                result = np.random.normal(0, temp, sample_number)[0]
        return result


    def accept_rate_scale(self, accept_list: list, scale=None) -> dict:
        ideal_rate = self.config.accept_rate()
        if scale is None:
            scale = {}
            scale['d_coor'] = 1.0
            scale['d_rot'] = 1.0
            scale['d_lat'] = 1.0
        # print(f"accept_list:{accept_list}")
        if len(accept_list) == 0:
            pass
        else:
            count_true = accept_list.count(True)
            accept_rate = round(count_true / len(accept_list), 3)
            if accept_rate > ideal_rate:
                for key in scale:
                    scale[key] *= 1.1
            elif accept_rate < ideal_rate:
                for key in scale:
                    scale[key] *= 0.9

        return scale
    

    def simulated_annealing(
        self,
        step: int | None = None,
        accept_list: list | None = None,
        count_number: None | int = None,
        best_energy: float | None = None,
        best_state: pyxtal | None = None,
        big_loop: int | None = None,
        d_lat: bool = False,
        random_sampling = 'normal',
        temp_final: float | None = None,
        temp_update: bool = True,
        accept_list_len: int = 50,
        write_output: bool = True,
        trajectory: bool = False,
        rate_scale: dict = None
    ):
        """
        !currently sampling scale has bit of issue and it's disabled.
        Inital parameters including:
        1. Current temperature
        2. Current structure (current state)
        3. Current step number
        4. Accept list to calculate acceptance rate
        5. Scale dict to adjust perturbation according to energy change with dof

        Set up perturbation to current state: 
            Random sampling * Scale dict * Scale dict 2 (according to acceptance rate)
        Compute objective function (energy)
        Determine to take new state or not 
        Update structure, temperature, step number and acceptance list.
        """
        # Initialize the starting state and best state
        max_step = self.config.internal_loop()
        temp_initial = self.config.temp_initial()
        cooling = self.config.temp_step()
        sub_command = self.config.sub_command()
        caltool = self.config.cal_tool()
        routine = self.config.routine_number()
        if routine == "compare":
            vasp_ml_compare = True
        else:
            vasp_ml_compare = False
        # Initialize best structure and best energy
        current_state = self.structure
        if best_energy is None:
            best_state = self.structure
            cal = CalculateEnergy(self.structure)
            best_energy = cal.energy_calculate(caltool, False, sub_command)[1]
        else:
            best_energy = best_energy
            best_state = best_state
        if accept_list is None:
            accept_list = []
        else:
            accept_list = accept_list
        # count for accept_list
        if count_number is not None:
            count = count_number
        else:
            count = 0
        # initialize step count
        if step is not None:
            step_initial = step + 1
            
        else:
            step_initial = 0
        
        # initialize temp
        if big_loop is not None:
            temp_current = self.cooling_rate(temp_initial, cooling, big_loop, temp_final)
        else:
            temp_current = temp_initial

        # scale_dict = self.sampling_scale(temp_current, current_state, d_lat=d_lat)
        
        for step in range(max_step):
            step_current = step_initial + step
            print(f"Step: {step_current}, Temperautre: {temp_current}")

            # Generate a neighboring state
            next_state = copy.deepcopy(current_state)
            next_state = apply_perturbation(next_state)
            # write trajectory
            if trajectory:
                trajectory_path = os.getcwd()+"/Trajectory"
                struc_py = structure_type_converter(next_state, "pymatgen")
                poscar = Poscar(struc_py)
                poscar.write_file(f"{trajectory_path}/POSCAR_{step}")
 
            # Compute the objective function value (energy) of the new state
            if vasp_ml_compare:
                cal_current = CalculateEnergy(current_state)
                cal_next = CalculateEnergy(next_state)
                current_energy_vasp = cal_current.energy_calculate("VASP", False, sub_command)[1]
                next_energy_vasp = cal_next.energy_calculate("VASP", False, sub_command)[1]
                next_energy_ml= cal_next.energy_calculate("CHGnet", False, sub_command)[1]

                current_energy = current_energy_vasp
                next_energy = next_energy_vasp

                print(f"Next state: VASP:{next_energy} eV, ML:{next_energy_ml} eV")
                # print(f"New state: {next_energy} eV")
                result_dir = os.getcwd() + "/Result"
                if not os.path.isfile(result_dir+"/Compare"):
                    with open(result_dir+"/Compare", 'w') as compare:
                        compare.write("Step, Temperautre/K, Energy-VASP/eV, Energy-ML/eV\n")
                with open(result_dir+"/Compare", 'a') as compare:
                    compare.write(f"{step_current}, {temp_current}, {next_energy}, {next_energy_ml}\n")
            
            else:
                cal_current = CalculateEnergy(current_state)
                cal_next = CalculateEnergy(next_state)
                print("----- Calculating energy of current state -----")
                current_cal = cal_current.energy_calculate(caltool, False, sub_command)
                current_energy = current_cal[1]
                current_energy_atom = current_cal[2]
                print("----- Calculating energy of new state -----")
                next_cal = cal_next.energy_calculate(caltool, False, sub_command)
                next_energy = next_cal[1]
                next_energy_atom = next_cal[2]
                print("-------------------------------------------------")
                print(f"Current state: {current_energy} eV, {current_energy_atom} eV/atom")
                print(f"New state: {next_energy} eV, {next_energy_atom} eV/atom")

            # Determine if we should accept the new state
            accept = False
            if next_energy < current_energy:
                current_state = next_state
                accept = True
                
            else:
                # Accept worse state with a probability depending on the temperature
                # accept_prob = 1.0 / (1.0 + math.exp(((next_energy - current_energy))/ temp_current))
                accept_prob = 1.0 / (1.0 + math.exp(((next_energy_atom - current_energy_atom))/ temp_current))
                if (accept_prob > random.uniform(0, 1) or accept_prob > 0.5) \
                    and next_energy_atom < 0 \
                    and next_energy/current_energy > 0.8:
                    current_state = next_state
                    accept = True
            print(f"Accept: {accept}")

            if len(accept_list) < accept_list_len:
                accept_list.append(accept)
            else:
                accept_list[count] = accept
                count += 1
                if count == accept_list_len:
                    count = 0
            count_true = accept_list.count(True)
            accept_rate = round(count_true / len(accept_list), 3)
            print(f"Acceptance rate: {round(accept_rate*100, 3)}%")
            print("\n")

            if next_energy < best_energy:
                best_state = next_state
                best_energy = next_energy

            # Write output
            if write_output:
                result_dir = os.getcwd() + "/Result"
                if not os.path.isfile(result_dir+"/Log"):
                    with open(result_dir+"/Log", 'w') as log:
                        log.write("Step, Temperautre/K, AcceptanceRate/%, CurrentEnergy/eV\n")
                with open(result_dir+"/Log", 'a') as log:
                    log.write(f"{step_current}, {temp_current}, {round(accept_rate*100., 2)}, {current_energy}\n")
                
                # if vasp_ml_compare:
                #     if not os.path.isfile(result_dir+"/Compare"):
                #         with open(result_dir+"/Compare", 'w') as compare:
                #             compare.write("Step, Temperautre, AcceptanceRate, Energy-VASP, Energy-ML\n")
                #     with open(result_dir+"/Compare", 'a') as compare:
                #         compare.write(f"{step_current}, {temp_current} K, {round(accept_rate*100., 2)}%, {current_energy_vasp} eV, {current_energy_ml} eV\n")

            # Updating temperature
            if temp_update:
                temp_current = self.cooling_rate(temp_initial, cooling, big_loop, temp_final)
            sys.stdout.flush()
        
        # count is for updating the accpet list
        return step_current, temp_current, accept_list, current_state, count, best_energy, best_state, rate_scale
