import os
import sys
import shutil
from readconfig import ReadConfig
from routine import Routine
from version import __version__

def print_logo():
    print(
        """
                  |\__/,|   (`
                _.|o o  |_   ) )
               -(((---(((--------
    Space Group Restricted Crystal Structure Prediction
        """
    )
    print("-------------------(version", __version__,")--------------------\n")
    print("A Python package for symmerty restricted crystal structure prediction")
    print("The source code is available at https://github.com/ColdSnaap/sgrcsp.git")
    print("Developed by ColdSnaap \n\n")
    sys.stdout.flush()


def clean():
    cal_dir = os.getcwd() + "/Calculation"
    if os.path.exists(cal_dir):
        shutil.rmtree(cal_dir)


if __name__ == "__main__":
    print_logo()
    config = ReadConfig()
    routine_number = config.routine_number()
    ready = config.ready()
    if ready:
        if routine_number == "routine1":
            optim = Routine()
            optim.routine1()
        elif routine_number == "routine2":
            optim = Routine()
            optim.routine2()
        elif routine_number == "routine3":
            optim = Routine()
            optim.routine3()
        elif routine_number == "compare":
            optim = Routine()
            optim.vasp_ml_compare()
        elif routine_number == "routine4":
            optim = Routine()
            optim.routine4()
        else:
            raise NameError(f"{routine_number} is not supported")
     