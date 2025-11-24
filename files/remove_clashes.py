import argparse
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm


class water:
    def __init__(self, OW, H1, H2):
        self.OW = OW
        self.H1 = H1
        self.H2 = H2


parser = argparse.ArgumentParser(
        prog='refine_results.py',
        description='a script that analyzes and refines dowser results')
parser.add_argument('-d', '--dowser_input', type=str,
                    default='PredictedInternal.pdb')
parser.add_argument('-s', '--solvent_input', type=str)
parser.add_argument('-o', '--output_pdb', type=str, default='no_clashes.pdb')


def record_unique_xyz(coords_info, i, coords_dict):
    coords_dict.__setitem__(coords_info, i)
    return coords_dict


def get_unique_dowser(coords_info, dowser_dict):
    coords_dict = {}
    for i in range(len(coords_info)):
        coords_dict = record_unique_xyz(coords_info[i], i, coords_dict)
    keys = list(coords_dict.keys())
    water_count = int(len(keys) / 3)
    for i in range(water_count):
        water_obj = water(keys[i * 3],
                          keys[i * 3 + 1], keys[i * 3 + 2])
        dowser_dict.__setitem__(water_obj, i + 1)

    dowser_unique = dict((v, k) for k, v in dowser_dict.items())
    return dowser_unique


def read_dowser_water(dowser_o):
    dowser_data_unique = {}
    dowser_file = dowser_o.readlines()
    dowser_file = [line for line in dowser_file if 'ATOM' or 'HETATM' in line]
    coords_info = [line[30:67] for line in dowser_file]
    print(f"There are {int(len(coords_info) / 3)} waters before removing duplicates...")
    dowser_unique = get_unique_dowser(coords_info, dowser_data_unique)
    unique_count = len(dowser_unique)
    print(f"There are {unique_count} waters after removing duplicates...")
    dowser_data = []
    keys = list(dowser_unique.keys())
    for i in range(len(keys)):
        current_water = dowser_unique[keys[i]]
        O_line = "ATOM  {:>5}".format(i % 99999 + 1) +\
                 "  OW  HOH A{:>4}    {}".format(1,
                                                 current_water.OW)
        H1_line = "ATOM  {:>5}".format(i % 99999 + 1) +\
                  "  H1  HOH A{:>4}    {}".format(1,
                                                  current_water.H1)
        H2_line = "ATOM  {:>5}".format(i % 99999 + 1) +\
                  "  H2  HOH A{:>4}    {}".format(1,
                                                  current_water.H2)

        one_line = O_line + H1_line + H2_line
        dowser_data.append(one_line)

    return dowser_data


def read_solvent_water(solvent_o):
    solvent_file = solvent_o.readlines()
    solvent_data = [line for line in solvent_file if 'SOL' and 'OW' in line]
    solvent_coords = np.array([np.array(line[30:56].split()).astype(float)
                               for line in solvent_data])

    return solvent_coords


def interaction_correction(dowser_E, n_neighbor, E_hbond=-2.5):
    if n_neighbor <= 4:
        dowser_E_corr = dowser_E + n_neighbor * E_hbond
    else:
        dowser_E_corr = dowser_E + 4 * E_hbond

    return dowser_E_corr


def check_hbond_neighbors(dowser_data, r=2.5, shell_thickness=0.5):
    neighbor_count = []
    dowser_E_corr = []
    water_num = np.array([int(x[7:12].strip()) for x in dowser_data])
    dowser_E = np.array([float(x[60:67]) for x in dowser_data])
    dowser_xyz = np.array([x[30:54].split() for
                           x in dowser_data]).astype(float)
    print("checking water neighbors...")
    for i in tqdm(range(len(dowser_data))):
        dist = np.sqrt(np.sum(np.square(dowser_xyz - dowser_xyz[i]), axis=1))
        within_shell = np.intersect1d(np.where(dist > r),
                                      np.where(dist <= (r + shell_thickness)))
        within_shell = within_shell.tolist()
        # print(f"water number {water_num[i]} has {len(within_shell)} neighbors")
        dowser_within_shell = [dowser_data[x] for x in within_shell]
        dowser_within_shell = sort_water_by_energy(dowser_within_shell)
        dowser_within_shell = remove_clashes(dowser_within_shell, r=2.5)
        n_neighbor = len(dowser_within_shell)
        # with open(f"neighbors/water_{water_num[i]}_neighbors.pdb", 'w') as neighbor_pdb:
        #     neighbor_pdb.writelines(dowser_within_shell)
        # print(f"water number {water_num[i]} has {n_neighbor} neighbors after removing clashes")
        neighbor_count.append(n_neighbor)
        dowser_E_corr.append(interaction_correction(dowser_E[i], n_neighbor, E_hbond=-2.5))
    neighbor_count = np.array(neighbor_count)
    dowser_E_corr = np.array(dowser_E_corr)
    water_info = pd.DataFrame({"water_num": water_num,
                               "Dowser_E": dowser_E,
                               "n_neighbor": neighbor_count,
                               "Dowser_E_corr": dowser_E_corr})
    # water_info.to_csv("water_info.csv", index=False)
    bins = np.arange(neighbor_count.min() - 0.5, neighbor_count.max() + 1.5, 1)
    plt.hist(neighbor_count, bins=bins)
    plt.xticks(np.arange(neighbor_count.min(), neighbor_count.max() + 1))
    plt.title("distribution of neighbor count within solvation shell")
    plt.xlabel("number of neighbors")
    plt.show()

    return water_info


def recheck_hbond_neighbors(dowser_data, water_info, r=2.5, shell_thickness=0.5):
    neighbor_count = []
    dowser_E_corr = []
    water_info["neighbors_after"] = "removed"
    water_num = np.array([int(x[7:12].strip()) for x in dowser_data])
    dowser_E = np.array([float(x[60:67]) for x in dowser_data])
    dowser_xyz = np.array([x[30:54].split() for
                           x in dowser_data]).astype(float)
    print("checking water neighbors...")
    for i in tqdm(range(len(dowser_data))):
        water_i = water_num[i]
        dist = np.sqrt(np.sum(np.square(dowser_xyz - dowser_xyz[i]), axis=1))
        within_shell = np.intersect1d(np.where(dist > r),
                                      np.where(dist <= (r + shell_thickness)))
        within_shell = within_shell.tolist()
        # print(f"water number {water_i} has {len(within_shell)} neighbors")
        dowser_within_shell = [dowser_data[x] for x in within_shell]
        dowser_within_shell = sort_water_by_energy(dowser_within_shell)
        dowser_within_shell = remove_clashes(dowser_within_shell, r=2.5)
        n_neighbor = len(dowser_within_shell)
        # with open(f"neighbors_after/water_{water_i}_neighbors.pdb", 'w') as neighbor_pdb:
        #     neighbor_pdb.writelines(dowser_within_shell)
        # print(f"water number {water_i} has {n_neighbor} neighbors after removing clashes")
        dowser_E_corr.append(interaction_correction(dowser_E[i], n_neighbor, E_hbond=-2.5))
        neighbor_count.append(n_neighbor)
        water_info.loc[water_info["water_num"] == water_i, "neighbors_after"] = n_neighbor
    neighbor_count = np.array(neighbor_count)
    dowser_E_corr = np.array(dowser_E_corr)
    water_info.to_csv("water_info_after.csv", index=False)
    bins = np.arange(neighbor_count.min() - 0.5, neighbor_count.max() + 1.5, 1)
    plt.hist(neighbor_count, bins=bins)
    plt.xticks(np.arange(neighbor_count.min(), neighbor_count.max() + 1))
    plt.title("distribution of neighbor count within solvation shell")
    plt.xlabel("number of neighbors")
    plt.show()

    return water_info


def remove_clashes(dowser_data, r: float = 2.5, E_threshold=-10):
    j = 1
    while True:
        # print(f"round {j}...")
        num_old = len(dowser_data)
        i = 0
        while i < len(dowser_data):
            dowser_xyz = np.array([x[30:54].split() for
                                   x in dowser_data]).astype(float)
            dowser_E = np.array([float(x[60:67]) for x in dowser_data])
            dist = np.sqrt(np.sum(np.square(dowser_xyz - dowser_xyz[i]), axis=1))
            clash_i = np.where(dist <= r)[0]
            clash_with_others = np.setdiff1d(clash_i, np.array(i))
            if clash_with_others.any():
                c = np.setdiff1d(clash_i, np.array(i))[0]
                # if dowser_E[c] >= E_threshold and dowser_E[i] >= E_threshold:
                # if the other water has lower energy
                if dowser_E[c] <= dowser_E[i]:
                    # remove current water
                    dowser_data.pop(i)
                else:
                    # remove the other water
                    dowser_data.pop(c)
            i += 1
        num_new = len(dowser_data)
        j += 1
        # print(f"{num_new} water molecules remain after checking clashes...")
        if num_old == num_new:
            break

    return dowser_data


def shuffle_water(dowser_data):
    dowser_data_randomized = random.sample(dowser_data, len(dowser_data))
    return dowser_data_randomized


def sort_water_by_energy(dowser_data):
    dowser_E = np.array([float(x[60:67]) for x in dowser_data])
    dowser_arr = np.array(dowser_data)
    dowser_data_sorted = dowser_arr[dowser_E.argsort()].tolist()

    return dowser_data_sorted


def add_mean_field_energy(dowser_data, dowser_E_corr):
    dowser_with_mean_field = []
    for i in range(len(dowser_data)):
        one_h2o = dowser_data[i]
        one_E = dowser_E_corr[i]
        dowser_E_new_str = f"{one_E:6.2f}"
        one_h2o_new = one_h2o[:60] + dowser_E_new_str + one_h2o[66:]

        dowser_with_mean_field.append(one_h2o_new)

    return dowser_with_mean_field


def energy_screening(dowser_data, E_cutoff=-4):
    total_E = 0
    dowser_within_cutoff = [x for x in dowser_data
                            if float(x[60:67]) < E_cutoff]
    for one_dowser in dowser_within_cutoff:
        total_E += float(one_dowser[60:67])

    return dowser_within_cutoff, total_E


if __name__ == "__main__":
    args = parser.parse_args()
    dowser_input = args.dowser_input
    solvent_input = args.solvent_input
    output_pdb = args.output_pdb

    with open(dowser_input, 'r') as dowser_o:
        dowser_data = read_dowser_water(dowser_o)
        # dowser_data = shuffle_water(dowser_data)
        dowser_data = sort_water_by_energy(dowser_data)

    with open(solvent_input, 'r') as solvent_o:
        solvent_coords = read_solvent_water(solvent_o)

    water_info = check_hbond_neighbors(dowser_data, r=2.5, shell_thickness=0.5)
    dowser_E_corr = water_info["Dowser_E_corr"]
    E_cutoff = -4  # kcal
    E_mean_field = -5  # kcal

    with open('./no_duplicates.pdb', 'w') as no_dupes:
        no_dupes.writelines(dowser_data)

    print(f"checking clashes of {len(dowser_data)} water molecules...")
    dowser_data = add_mean_field_energy(dowser_data, dowser_E_corr)
    no_clash = remove_clashes(dowser_data, r=2.5, E_threshold=-10)
    num_new = len(no_clash)
    # no_clash, total_E = energy_screening(no_clash, E_cutoff=E_cutoff)
    print(f"{num_new} water molecules remain after energy screening...")
    water_info = recheck_hbond_neighbors(no_clash, water_info=water_info, r=2.5,
                                         shell_thickness=0.5)
    # print(f"total energy is {total_E:.2f} kcal")
    with open(output_pdb, 'w') as output:
        output.writelines(no_clash)
    pass
