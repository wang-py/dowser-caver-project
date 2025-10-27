import argparse
import numpy as np


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
        O_line = "ATOM  {:>5}".format(i + 1) +\
                 "  OW  HOH A{:>4}    {}".format(i + 1,
                                                 current_water.OW)
        H1_line = "ATOM  {:>5}".format(i + 1) +\
                  "  H1  HOH A{:>4}    {}".format(i + 1,
                                                  current_water.H1)
        H2_line = "ATOM  {:>5}".format(i + 1) +\
                  "  H2  HOH A{:>4}    {}".format(i + 1,
                                                  current_water.H2)

        one_line = O_line + H1_line + H2_line
        dowser_data.append(one_line)

    return dowser_data


def remove_clashes(dowser_data, r: float = 2.5, E_threshold=-10):
    i = 0
    while i < len(dowser_data):
        dowser_xyz = np.array([x[30:54].split() for
                               x in dowser_data]).astype(float)
        dowser_E = np.array([float(x[60:66]) for x in dowser_data])
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

    return dowser_data


def add_mean_field_energy(dowser_data, shell_radius=3, E_mean_field=-5):
    dowser_xyz = np.array([x[30:54].split() for
                           x in dowser_data]).astype(float)
    dowser_with_mean_field = []
    for i in range(len(dowser_data)):
        one_h2o = dowser_data[i]
        one_xyz = np.array(one_h2o[30:54].split()).astype(float)
        one_E = float(one_h2o[60:66].strip())
        dist = np.sqrt(np.sum(np.square(dowser_xyz - one_xyz), axis=1))
        within_shell = np.where(dist <= shell_radius)[0]
        h2o_within_shell = np.setdiff1d(within_shell, np.array(i))

        if h2o_within_shell.any():
            dowser_E_new = one_E + E_mean_field
            dowser_E_new_str = f"{dowser_E_new:6.2f}"
            one_h2o_new = one_h2o[:60] + dowser_E_new_str + one_h2o[66:]
        else:
            one_h2o_new = one_h2o

        dowser_with_mean_field.append(one_h2o_new)

    return dowser_with_mean_field


def energy_screening(dowser_data, E_cutoff=-4):
    dowser_within_cutoff = [x for x in dowser_data
                            if float(x[60:66]) < E_cutoff]

    return dowser_within_cutoff


if __name__ == "__main__":
    args = parser.parse_args()
    dowser_input = args.dowser_input
    output_pdb = args.output_pdb

    with open(dowser_input, 'r') as dowser_o:
        dowser_data = read_dowser_water(dowser_o)

    E_cutoff = -4  # kcal
    E_mean_field = -5  # kcal

    with open('./no_duplicates.pdb', 'w') as no_dupes:
        no_dupes.writelines(dowser_data)

    print(f"checking clashes of {len(dowser_data)} water molecules...")
    no_clash = dowser_data
    i = 1
    while True:
        print(f"round {i}...")
        num_old = len(no_clash)
        no_clash = remove_clashes(no_clash, r=2.5, E_threshold=-10)
        num_new = len(no_clash)
        i += 1
        print(f"{num_new} water molecules remain after checking clashes...")
        if num_old == num_new:
            break
    no_clash = add_mean_field_energy(no_clash, shell_radius=3,
                                     E_mean_field=E_mean_field)
    no_clash = energy_screening(no_clash, E_cutoff=E_cutoff)
    print(f"{num_new} water molecules remain after energy screening...")
    with open(output_pdb, 'w') as output:
        output.writelines(no_clash)
    pass
