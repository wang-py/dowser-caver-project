import argparse

parser = argparse.ArgumentParser(
        prog='refine_results.py',
        description='a script that analyzes and refines dowser results')
parser.add_argument('-d', '--dowser_input', type=str,
                    default='PredictedInternal.pdb')
parser.add_argument('-o', '--output_pdb', type=str, default='unique.pdb')


def read_dowser_water(dowser_o):
    dowser_file = dowser_o.readlines()
    dowser_data_unique = list(set([line[30:66] for line in dowser_file
                                   if 'ATOM' and ' O ' in line]))
    dowser_data = []
    for i in range(len(dowser_data_unique)):
        # one_line = "ATOM" + "  " + f"{i + 1}".rjust(4) + "  O  " + " HOH " +\
        #     "A " + f" {i + 1}  ".rjust(3) + "   " + dowser_data_unique[i] +\
        #     "\n"
        one_line = "ATOM  {:>5}".format(i + 1) +\
                "  O   HOH A{:>4}    {}\n".format(i + 1,
                                                  dowser_data_unique[i])
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
