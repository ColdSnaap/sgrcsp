import ast
import os
from molecule import MolecularObject

class ReadConfig:
    """
    Read configuration from a file and apply default values for any missing parameters.
    """
    def __init__(self):
        self.config = {}
        self.input_dir = os.getcwd() + "/Input"
        self.input_file_path = self.input_dir + "/input.txt"
        # Define default values
        default_values = {
        "SpaceGroup": 1,
        "ForceMax": 0.1,
        "ScalePerturbate": False,
        "Optimization": "routine1",
        "TempStep": "fast",
        "TempInital": 2.0,
        "TempFinal": 0.1,
        "Ready": False,
        "SubCommand": None,
        "MinDistance": 1.4,
        "Volume": 0.6,
        "MaxComb": 1000,
        "MaxCombCal": 20,
        "InternalLoop": 400,
        "AcceptRate": 0.25,
        "StepMax": 80000,
        "VaspMlCompare": False,
        "InitialRelax": False,
        "Continue": False
        }
        # Try to open and read the file
        try:
            with open(self.input_file_path, 'r') as file:
                for line in file:
                    # Ignore empty lines and lines starting with '#'
                    if not line or line.startswith('#'):
                        continue
                    # Ignore empty lines and lines without '='
                    if "=" not in line:
                        continue
                    key, value = line.strip().split('=', 1)
                    self.config[key.strip()] = value.strip()
        except FileNotFoundError:
            print("Configuration file not found.")
            return None
        except Exception as e:
            print(f"An error occurred: {e}")
            return None
        
        # Apply default values for any missing parameters
        for key, value in default_values.items():
            if key not in self.config:
                self.config[key] = value

        
    def mol_list(self):     
        if "MolFiles" in self.config:
            mol_list = []
            mol_number = len(self.config["MolFiles"].split())
            for mol_index in range(mol_number):
                mol_iterate = self.config["MolFiles"].split()[mol_index]
                mol = MolecularObject(self.input_dir + f"/{mol_iterate}").pymatgen()
                mol_list.append(mol)
        return mol_list
    
    
    def space_group(self):
        if "SpaceGroup" in self.config and not isinstance(self.config["SpaceGroup"], int):
            sg_list = []
            for item in self.config["SpaceGroup"].split():
                if '-' in item:  
                    start, end = item.split('-')  
                    sg_list.append([int(start), int(end)]) 
                else:
                    sg_list.append([int(item)])
        else:
            sg_list = [[self.config["SpaceGroup"]]]
        sg_list_int = []
        for sg_item in sg_list:
            if len(sg_item) == 1:
                sg_list_int.append(int(sg_item[0]))
            else:
                for i in range(int(sg_item[0]), int(sg_item[1])+1):
                    sg_list_int.append(i)
        sg_list_int.sort()
        unique_lst = list(dict.fromkeys(sg_list_int))
        return unique_lst
    

    def mol_number(self):
        if "MolNumber" in self.config:
            mol_numner = []
            mol_number_str = self.config["MolNumber"].split()
            mol_number = [int(item) for item in mol_number_str]
        return mol_number
    

    def max_multiplier(self):
        max_muti = int(self.config["RatioMaxMultiplier"])
        return max_muti

    
    def mol_sites(self):
        mol_sites = ast.literal_eval(self.config["MolSites"])
        return mol_sites


    def cal_tool(self):
        cal_tool = self.config["Caltool"]
        return cal_tool
    
    
    def fmax(self):
        fmax = self.config["ForceMax"]
        return float(fmax)
    

    def sub_command(self):
        sub_command = self.config["SubCommand"]
        return str(sub_command)
    
    
    def temp_initial(self):
        temp = self.config["TempInital"]
        return float(temp)
    

    def internal_loop(self):
        loop_step = self.config["InternalLoop"]
        return int(loop_step)

    
    def max_step(self):
        max_step = self.config["StepMax"]
        return int(max_step)


    def accept_rate(self):
        accept_rate = self.config["AcceptRate"]
        return float(accept_rate)


    def routine_number(self):
        routine = self.config["Optimization"]
        return routine
    

    def max_comb(self):
        """
        The biggest number allowed for combination
        """
        max_comb = self.config["MaxComb"]
        return int(max_comb)
    

    def ready(self):
        ready = self.config["Ready"]
        actual_boolean = ready == "True"
        return actual_boolean
    

    def max_com_cal(self):
        max_com_cal = self.config["MaxCombCal"]
        return int(max_com_cal)
    

    def temp_step(self):
        temp_step = self.config["TempStep"]
        return temp_step

    
    def job_continue(self):
        jc = self.config["Continue"]
        actual_boolean = jc == "True"
        return actual_boolean
    

    def min_dis(self):
        min_dis = self.config["MinDistance"]
        return float(min_dis)

    
    def volume(self):
        volume = self.config["Volume"]
        return float(volume)

    
    def compare(self):
        compare = self.config["VaspMlCompare"]
        actual_boolean = compare == "True"
        return actual_boolean

    
    def initial_relax(self):
        relax = self.config["InitialRelax"]
        actual_boolean = relax == "True"
        return actual_boolean