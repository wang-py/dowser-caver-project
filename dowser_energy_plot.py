import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
import sys
import re

def get_energy_from_xvg(input_xvg):
    with open(input_xvg, 'r') as xvg:
        data_str = [line.split()[1:] for line in xvg if '#' not in line and '@' not in line]
    
    data = np.array(data_str).astype(float)

    return data

def get_site_number_from_energy_file(energy_file):
    filename = os.path.basename(energy_file)
    sn = re.findall('\_\d+\_', filename)
    sn = str(sn[0])
    sn = sn.strip('_')
    site_number = float(sn)
    #site_number = float(filename.split('_')[5])

    return site_number
    
def get_energy_vs_site(energy_files):
    """
    function that gets individual average energies from input xvgs
    """
    num_of_pts = len(energy_files)
    total_energies = np.zeros(num_of_pts)
    site_number = np.zeros(num_of_pts)
    for i in range(num_of_pts):
        energy = get_energy_from_xvg(energy_files[i])
        site_number[i] = get_site_number_from_energy_file(energy_files[i])
        coulomb = energy[:, [0,2]]
        LJ = energy[:, [1,3]]
        avg_total_energy = np.mean(coulomb + LJ)
        total_energies[i] = avg_total_energy

    return total_energies, site_number

def plot_energy_vs_site(dowser_1, dowser_2, sites, output_filename):
    cal_to_joules = 4.184
    label_fontsize=16
    fig, ax = plt.subplots(figsize=(14,7))
    ax.set_xticks(sites)
    ax.tick_params(axis='x', labelsize=label_fontsize)
    ax.tick_params(axis='y', labelsize=label_fontsize)
    ax.plot(sites, dowser_1, 'rv', label='dowser refinement run', markersize=10)
    ax.plot(sites, dowser_2, 'bv', label='dowser refinement run with EM', markersize=10)
    #plt.title("total energy vs site number", fontsize=label_fontsize)
    bulk_energy = -42 / cal_to_joules
    plt.axhline(bulk_energy, color='k', linestyle='--', label='energy of water in bulk %.1f kCal/mol'%bulk_energy)
    plt.xlabel("site number", fontsize=label_fontsize)
    plt.ylabel("energy [kCal/mol]", fontsize=label_fontsize)
    plt.legend(loc="best")
    plt.savefig(output_filename+".png", dpi=200)
    plt.show()
    pass

def plot_gmx_dowser_energy_vs_site(total_energies, sites, dowser_energies, dowser_hit_stats):
    cal_to_joules = 4.1868
    fig, ax = plt.subplots()
    plt.plot(sites, total_energies / cal_to_joules, 'bo', label='gromacs')
    ax.set_xticks(sites)
    dowser_hits = sites[dowser_hit_stats].astype(int)
    overpredictions = sites[~dowser_hit_stats].astype(int)
    plt.scatter(dowser_hits, dowser_energies[dowser_hits-1], marker='o', label='dowser hits', color='green')
    plt.scatter(overpredictions, dowser_energies[overpredictions-1], marker='o', label='dowser overpredictions', color='orange')
    plt.title("total energy vs site number")
    bulk_energy = -42 / cal_to_joules
    plt.axhline(bulk_energy, color='k', linestyle='--', label='energy of water in bulk %.1f kCal/mol'%bulk_energy)
    plt.xlabel("site number")
    plt.ylabel("energy (Coulomb + LJ) [kCal/mol]")
    plt.legend()
    plt.show()
    pass

def get_dowser_energies(dowser_energy_file):
    with open(dowser_energy_file, 'r') as DE:
        dowser_energies = [float(line) for line in DE.readlines()]
    return np.array(dowser_energies)

if __name__ == "__main__":
    dowser_1 = sys.argv[1]
    dowser_2 = sys.argv[2]
    if len(sys.argv) > 3:
        output_fig_filename = sys.argv[3]
    else:
        output_fig_filename = "output_fig"
    dowser_energies_1 = get_dowser_energies(dowser_1)
    dowser_energies_2 = get_dowser_energies(dowser_2)
    sites = np.arange(1, len(dowser_energies_1)+1)
    plot_energy_vs_site(dowser_energies_1, dowser_energies_2, sites, output_fig_filename)

    pass