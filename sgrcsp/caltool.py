import os
import subprocess
import shutil
import sys
from pathlib import Path
from pymatgen.core import Structure
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import StructOptimizer
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp import Poscar
from structure import structure_type_converter

class CalculateEnergy:
    
    def __init__(self, structure, molecule=True) -> None:
        
        self.structure = structure_type_converter(structure, "pymatgen", molecule)


    def vasp_input(self, relax: bool):
        input_dir = os.getcwd() + "/Input"
        vasp_relax_dir = os.getcwd() + "/StrucRelax"
        vasp_en_dir = os.getcwd() + "/EnergyCal"
        poscar = Poscar(self.structure)

        # generalte POTCAR
        potcar_file = []
        element_list = self.structure.elements
        for element in element_list:
            potcar_ele = input_dir+f"/POTCAR_{str(element)}"
            potcar_file.append(potcar_ele)
        
        with open(input_dir+"/POTCAR", 'w') as outfile:
            for file_path in potcar_file:
                with open(file_path, 'r') as infile:
                    # Read and write line by line
                    for line in infile:
                        outfile.write(line)

        if relax:
            poscar.write_file(vasp_relax_dir+'/POSCAR')
            shutil.copy(input_dir+"/INCAR_relax", vasp_relax_dir+"/INCAR")
            shutil.copy(input_dir+"/KPOINTS_relax", vasp_relax_dir+"/KPOINTS")
            shutil.copy(input_dir+"/POTCAR", vasp_relax_dir+"/POTCAR")
        else:
            poscar.write_file(vasp_en_dir+'/POSCAR')
            shutil.copy(input_dir+"/INCAR", vasp_en_dir+"/INCAR")
            shutil.copy(input_dir+"/KPOINTS", vasp_en_dir+"/KPOINTS")
            shutil.copy(input_dir+"/POTCAR", vasp_en_dir+"/POTCAR")


    def energy_calculate(
        self, 
        tool: str,
        relax: bool,
        sub_command: str = None,
        fmax: float = 0.01
    ):

        if tool == 'VASP':
            if relax:
                self.vasp_input(relax=True)
            else:
                self.vasp_input(relax=False)
            structure, energy, energy_atom = self.vasp(sub_command, relax)
        elif tool == 'CHGnet':
            structure, energy, energy_atom = self.mlp_CHGNet(fmax, relax)
        
        return structure, energy, energy_atom    


    def vasp(self, submit_command: str, relax: bool = False) -> list:
        """
        If relax = False, relax_structure = False 
        """
        # root directory 
        root = os.getcwd()
        # determine working dir
        if relax:
            dir = root+"/StrucRelax"
        else:
            dir = root+"/EnergyCal"
        # change working directory
        # poscar = Poscar(self.structure)
        # poscar.write_file(dir+'/POSCAR')
        os.chdir(dir)
        structure_final = self.structure
        
        print("VASP running ...", flush=True)

        vasp_run = subprocess.run(submit_command, shell=True, capture_output=True, text=True)
        if vasp_run.returncode == 0:
            print(f"VASP(relax={relax}) job has completed successfully.", flush=True)
        #     print("Output:", vasp_run.stdout)
        else:
            print("Error in VASP job:", vasp_run.stderr)
        outcar_path = os.getcwd() + "/OUTCAR"
        outcar = Outcar(outcar_path)
        energy = round(outcar.final_energy, 5)
        atom_number = len(self.structure)
        if relax:
            structure_contcar = Poscar.from_file(os.getcwd() + "/CONTCAR")
            structure_final = structure_contcar.structure
        energy_per_atom = round(outcar.final_energy / atom_number, 5)

        os.chdir(root)

        return structure_final, energy, energy_per_atom


    def mlp_CHGNet(self, fmax: float = 0.1, relax: bool = False) -> list:
        relax_structure = False
        chgnet = CHGNet.load()
        if relax==False:
            # print(self.structure)
            prediction = chgnet.predict_structure(self.structure)
            energy = round(prediction['e'] * len(self.structure), 5)
            energy_per_atom = round(float(prediction['e']), 5)
        elif relax==True:
            relaxer = StructOptimizer()
            result = relaxer.relax(self.structure, fmax=fmax)
            energy = round(result['trajectory'].energies[-1], 5)
            relax_structure = Poscar(result['final_structure'])
            energy_per_atom = round(result['trajectory'].energies[-1] / len(self.structure), 5)

        return relax_structure, energy, energy_per_atom