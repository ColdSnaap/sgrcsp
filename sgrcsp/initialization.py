import os
import re
import shutil
from readconfig import ReadConfig
from combination import WyckoffCombinations

class Initializer:
    
    def __init__(self, space_group: int) -> None:
        self.config = ReadConfig()
        self.mol_list = self.config.mol_list()
        self.mol_number = self.config.mol_number()
        self.max_multiplier = self.config.max_multiplier()
        self.max_comb = self.config.max_comb()
        self.max_comb_cal = self.config.max_com_cal()
        self.space_group = space_group

        self.cal_dir = os.getcwd()+'/Calculation'
        os.makedirs(self.cal_dir, exist_ok=True)
    

    def write_comb_info(self, comb: list, file_path: str = os.getcwd()) -> None:
        if type(comb[0][0]) is str:
            with open(file_path+'/input.txt', 'a') as file:
                file.write(f'\nMolSites = {comb}')
        else:
            with open(file_path+'/combination_info', 'w') as file:
                for i, item in enumerate(comb):
                    file.write(f'{i+1}: {item}\n')


    def write_ready(self, file_path) -> None:
        with open(file_path, 'a') as file:
            file.write(f'\nReady = True')

    
    def modify_update(self, file_path, key, new_value):
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
    

    def mol_number_from_list(self, lists):
        result = []
        # Iterate over each sublist in the main list
        for sublist in lists:
            total = 0
            # Iterate over each string in the sublist
            for item in sublist:
                # Use regular expressions to find all digits up to the first non-digit character
                number = re.match(r'\d+', item)
                if number:
                    total += int(number.group())  # Convert the number to int and add to the total
            # Append the sum of this sublist to the result list
            result.append(total)
        numbers_str = " ".join(str(num) for num in result)
        return numbers_str


    def make_comb_dir(self, ratio: bool = True) -> None:
        comb = WyckoffCombinations(self.mol_list)
        if ratio:
            comb_list = comb.mol_ratio_comb_list_sg(self.mol_number, self.space_group, self.max_multiplier)
        else:
            comb_list = comb.mol_numb_comb_list_sg(self.mol_number, self.space_group, max_comb=self.max_comb)
        if len(comb_list[0]) != 0 and len(comb_list[0]) <= self.max_comb_cal:
            group_dir = f'{self.cal_dir}/group{self.space_group}'
            if os.path.exists(group_dir):
                shutil.rmtree(group_dir)
            os.makedirs(group_dir)
            self.write_comb_info(comb_list[0], group_dir)
            for i, item in enumerate(comb_list[0]):
                case_dir = f'{group_dir}/{i+1}'
                if os.path.exists(case_dir):
                    shutil.rmtree(case_dir)
                os.makedirs(case_dir)
                os.makedirs(case_dir + '/StrucRelax')
                os.makedirs(case_dir + '/EnergyCal')
                os.makedirs(case_dir + '/Result')
                shutil.copytree(os.getcwd()+'/Input', case_dir+'/Input')

                for filename in os.listdir(os.getcwd()):
                    if filename.endswith('.py'):
                        source_file = os.path.join(os.getcwd(), filename)
                        destination_file = os.path.join(case_dir, filename)

                        # Copy the file to the destination directory
                        shutil.copy2(source_file, destination_file)

                self.write_comb_info(item, case_dir+'/Input')
                self.modify_update(case_dir+'/Input/input.txt', "SpaceGroup", self.space_group)
                # Get mol_number
                mol_number = self.mol_number_from_list(item)
                self.modify_update(case_dir+'/Input/input.txt', "MolNumber", mol_number)
                self.write_ready(case_dir+'/Input/input.txt')

            # print(f"Space Group ")