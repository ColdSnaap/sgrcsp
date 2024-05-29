import sys
import os
from pyxtal.symmetry import Group
from readconfig import ReadConfig
from pyxtal.molecular_crystal import molecular_crystal

class WyckoffCombinations:
    
    def __init__(self, mol_list: list) -> None:
        config = ReadConfig()
        self.mol_list = mol_list
        self.max_comb = config.max_comb()


    def mol_numb_comb_list_sg(self, mol_number: list, sg: int) -> list:
        """
        wyckoff combinations with specific molecule number in a single space group
        """
        result_list = []
        pyxtal_sg = Group(sg)
        comb = pyxtal_sg.list_wyckoff_combinations(mol_number)[0]
        for i, wyckoff_list in enumerate(comb):
            mol_crystal = molecular_crystal(dim=3, group=sg, molecules=self.mol_list, numMols=mol_number, factor=200.0, sites=wyckoff_list)
            if mol_crystal.valid:
                result_list.append(wyckoff_list)
        comb_number = len(result_list)
        return result_list, comb_number

    
    def mol_ratio_comb_list_sg(self, mol_ratio: list, sg: int, max_multiplier) -> list:
        """
        wyckoff combinations with specific molecule ratio in a single space group
        """
        result_list = []
        pyxtal_sg = Group(sg)
        symbol = pyxtal_sg.symbol
        hit_max = False

        for muti in range(1, max_multiplier + 1):
            mol_number = [int(item) * muti for item in mol_ratio]
            comb = pyxtal_sg.list_wyckoff_combinations(mol_number)[0]
            if len(comb) != 0:
                valid_number = 0
                for i, wyckoff_list in enumerate(comb):
                    try:
                        mol_crystal = molecular_crystal(dim=3, group=sg, molecules=self.mol_list, numMols=mol_number, factor=200.0, sites=wyckoff_list) 
                    except IndexError:
                        pass
                    if mol_crystal.valid:
                        result_list.append(wyckoff_list)
                        valid_number += 1
                    
                    sys.stdout.write(f"\r{symbol} ({sg}), Z = {muti}, checking: {i+1}/{len(comb)}, valid: {valid_number}, sites: {wyckoff_list}                      ")
                    sys.stdout.flush()

                    if len(result_list) >= self.max_comb:
                        hit_max = True
                        break 
                if len(result_list) >= self.max_comb:
                    hit_max = True
                    break
                sys.stdout.write("\n")
        comb_number = len(result_list)

        # write the output file
        comb_result = os.getcwd()+"/comb.txt"
        if not os.path.exists(comb_result):
            with open("comb.txt", 'w') as file:
                file.write("sg,comb\n")
        with open(comb_result, 'a') as comb:
            if hit_max:
                comb.write(f"{sg},>{self.max_comb}\n")
            elif comb_number != 0:
                comb.write(f"{sg},{comb_number}\n")
        
        return result_list, comb_number
