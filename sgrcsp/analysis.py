import os
import shutil
import matplotlib.pyplot as plt
import pandas as pd
from optimization import GlobalOptimize


def vasp_ml_compare(csv_file, title=None, write_file=False, data_number=100):
    """
    retrun accuarcy of universal potential
    """
    # Read the data from the CSV file
    # df = pd.read_csv('/Users/qizhang/Desktop/paper_new/VASP_ML_compare/PSLi_compare.csv')
    df = pd.read_csv(csv_file)
    df = df[df['Step'] <= data_number]
    # check accuracy 
    df['Change_Sign'] = 0
    df['Difference'] = 0
    # Perform the subtraction and determine the sign change
    for i in range(1, len(df)):
        previous_vasp = df.at[i-1, ' Energy-VASP/eV']
        current_vasp = df.at[i, ' Energy-VASP/eV']
        previous_ml = df.at[i-1, ' Energy-ML/eV']
        current_ml = df.at[i, ' Energy-ML/eV']
        
        vasp_diff = current_vasp - previous_vasp
        ml_diff = current_ml - previous_ml

        # Check if the sign of the difference has changed
        if (vasp_diff * ml_diff) < 0:
            df.at[i, 'Change_Sign'] = 1
        else:
            dif =  abs(vasp_diff-ml_diff)/abs(vasp_diff)
            # print(f"{vasp_diff}, {ml_diff}")
            # print(dif)
            df.at[i, 'Difference'] = dif
    # print(df)
    count_change_sign_1 = df['Change_Sign'].sum()
    
    count = 0
    average_sum = 0
    for i in range(len(df)):
        if df.at[i, 'Difference'] != 0:
            average_sum += df.at[i, 'Difference']
            count += 1
    
    if count == 0:
        raise ValueError("count is zero")
    average = average_sum / count
    # print(df.columns)
    # Plot the data
    plt.figure(figsize=(7, 5))
    plt.plot(df['Step'], df[' Energy-VASP/eV'], label=' VASP', marker='o')
    plt.plot(df['Step'], df[' Energy-ML/eV'], label=' CHGnet', marker='x')

    plt.xlabel('Step', fontsize=14)
    plt.ylabel('Energy (eV)', fontsize=14)
    # plt.title('Step vs Energy ($Li_3PS_4$)', fontsize=16)
    plt.title(title, fontsize=16)
    plt.legend()
    plt.grid(False)
    # plt.show()

    if write_file:
        plt.savefig('fig.eps', format='eps')

    return plt, df, count_change_sign_1, average


def check_trajectory(structure, perturb_dict):
    trajectory_path = os.getcwd()+"/Trajectory"
    if os.path.isdir(trajectory_path):
        shutil.rmtree(trajectory_path)
    os.mkdir(trajectory_path)
    op = GlobalOptimize(structure)
    op.simulated_annealing(trajectory=True, perturb_dict=perturb_dict, write_output=False, temp_update=False)