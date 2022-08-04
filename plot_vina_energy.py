import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
from glob import glob

def get_energy_of_one_prediction(pdbqt):
    with open(pdbqt,'r') as f:
        data = [line for line in f if "REMARK INTER + INTRA" in line]
    
    num_of_predictions = len(data)
    energies = np.zeros((num_of_predictions, 1))
    for i in range(num_of_predictions):
        line  = data[i]
        energies[i] = float(line.split()[4])
    
    return energies

def plot_all_energies(all_energies):
    prediction_index = np.array([])
    counter = 0
    for energy in all_energies:
        energy_len = len(energy)
        prediction_index = np.append(prediction_index, np.repeat(counter, energy_len))
        counter += 1
    all_energies = np.concatenate(all_energies)
    plt.figure()
    plt.title("AutoDock Vina Energies Per Site")
    plt.ylabel("Energy [kcal/mol]")
    plt.ylim([-3, 3])
    plt.xlabel("Docking site number")
    plt.scatter(prediction_index, all_energies)
    plt.show()

if __name__ == "__main__":
    prediction_path = sys.argv[1]
    prediction_files = sorted(glob(prediction_path + "/point_*.pdbqt"), key=os.path.getmtime)
    all_energies = []
    for pdbqt in prediction_files:
        one_energy = get_energy_of_one_prediction(pdbqt)
        all_energies.append(one_energy)

    plot_all_energies(all_energies)
