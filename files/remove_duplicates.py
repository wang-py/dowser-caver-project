import argparse


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
parser.add_argument('-o', '--output_pdb', type=str, default='unique.pdb')


def get_unique_dowser(coords_info, dowser_dict):
    water_count = int(len(coords_info) / 3)
    for i in range(water_count):
        water_obj = water(coords_info[i * 3],
                          coords_info[i * 3 + 1], coords_info[i * 3 + 2])
        dowser_dict.__setitem__(water_obj, i + 1)

    dowser_unique = dict((v, k) for k, v in dowser_dict.items())
    return dowser_unique


def read_dowser_water(dowser_o):
    dowser_data_unique = {}
    dowser_file = dowser_o.readlines()
    dowser_file = [line for line in dowser_file if 'ATOM' in line]
    coords_info = [line[30:67] for line in dowser_file]
    print(f"There are {int(len(coords_info) / 3)} waters before removing duplicates...")
    dowser_unique = get_unique_dowser(coords_info, dowser_data_unique)
    unique_count = int(len(dowser_data_unique) / 3)
    print(f"There are {unique_count} waters after removing duplicates...")
    dowser_data = []
    keys = list(dowser_unique.keys())
    for i in range(len(keys)):
        current_water = dowser_unique[keys[i]]
        O_line = "ATOM  {:>5}".format(i + 1) +\
                 "  O   HOH A{:>4}    {}".format(i + 1,
                                                 current_water.OW)
        H1_line = "ATOM  {:>5}".format(i + 1) +\
                  "  H   HOH A{:>4}    {}".format(i + 1,
                                                  current_water.H1)
        H2_line = "ATOM  {:>5}".format(i + 1) +\
                  "  H   HOH A{:>4}    {}".format(i + 1,
                                                  current_water.H2)

        one_line = O_line + H1_line + H2_line
        dowser_data.append(one_line)

    return dowser_data


if __name__ == "__main__":
    args = parser.parse_args()
    dowser_input = args.dowser_input
    output_pdb = args.output_pdb

    with open(dowser_input, 'r') as dowser_o:
        dowser_data = read_dowser_water(dowser_o)
    with open(output_pdb, 'w') as output:
        output.writelines(dowser_data)
    pass
