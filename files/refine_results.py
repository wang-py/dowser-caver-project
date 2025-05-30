import numpy as np
import os
import subprocess
import copy
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(
        prog='refine_results.py',
        description='a script that analyzes and refines dowser results')
parser.add_argument('-d', '--dowser_input', type=str)
parser.add_argument('-s', '--structure_input', type=str)
parser.add_argument('-o', '--refined_pdb', type=str, default='refined.pdb')
parser.add_argument('-c', '--cutoff', type=float, default=-4.0)


def water_hetatm_replacement(position, structure_data):
    structure_data[position] = structure_data[position].replace('ATOM  ', 'HETATM')
    structure_data[position+1] = structure_data[position+1].replace('ATOM  ', 'HETATM')
    structure_data[position+2] = structure_data[position+2].replace('ATOM  ', 'HETATM')
    return structure_data


# ./reform -pdbin lipolytica_JKAHN_charged_with_dowser_waters_HOH23_HETATM.pdb
# -pdbout lipolytica_JKAHN_charged_with_dowser_waters_HOH23_HETATM_DOWSER.pdb
def run_reform(pdbin, pdbout):
    reform_args = ("./reform", "-pdbin", pdbin, "-pdbout", pdbout)
    popen = subprocess.Popen(reform_args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    popen.wait()


# ./placeWat lipolytica_JKAHN_charged_with_dowser_waters_HOH23_HETATM_DOWSER.pdb
# lipolytica_JKAHN_O23.pdb rotate > lipolytica_JKAHN_charged_with_dowser_waters
# HOH23_HETATM_docking_1.pdb
def run_placeWat(dowser_pdb, current_water_pdb, output):
    placeWat_args = ('./placeWat', dowser_pdb, current_water_pdb, 'rotate')
    # print(placeWat_args)
    popen = subprocess.Popen(placeWat_args, stdout=output, stderr=subprocess.DEVNULL)
    popen.wait()


def remove_disqualified_water(re_eval_pdb, cutoff):
    with open(re_eval_pdb, 'r') as rep:
        re_eval_results = rep.readlines()
    num_of_atoms = len(re_eval_results)
    refined = []
    for i in tqdm(range(0, num_of_atoms, 3)):
        # print(f"checking water {int(i / 3 + 1)}...", end='\n\r')
        energy_after_EM = float(re_eval_results[i][60:66])
        if energy_after_EM < cutoff:
            refined.append(re_eval_results[i])
            refined.append(re_eval_results[i+1])
            refined.append(re_eval_results[i+2])
            print(f"water {int(i / 3 + 1)} has an energy of" +
                  f" {energy_after_EM:.2f} kCal < {cutoff:.2f} kCal," +
                  " thus it was kept.")
        else:
            print(f"water {int(i / 3 + 1)} has an energy of" +
                  f" {energy_after_EM:.2f} kCal > {cutoff:.2f} kCal," +
                  " thus it was removed.")
    return refined


def read_dowser_water(dowser_o):
    dowser_file = dowser_o.readlines()
    dowser_data_unique = list(set([line[30:66] for line in dowser_file
                                   if 'ATOM' and ' O ' in line]))
    dowser_data = []
    for i in range(len(dowser_data_unique)):
        one_line = "ATOM" + "    " + f" {i + 1} ".rjust(4) + " O " + " HOH " +\
            "A " + f" {i + 1}  ".rjust(3) + "   " + dowser_data_unique[i] +\
            "         " + " O \n"
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
            if dowser_E[c] >= E_threshold and dowser_E[i] >= E_threshold:
                # if the other water has lower energy
                if dowser_E[c] <= dowser_E[i]:
                    # remove current water
                    dowser_data.pop(i)
                else:
                    # remove the other water
                    dowser_data.pop(c)
        i += 1

    return dowser_data


def update_structure_data(dowser_data, structure_data):
    protein_data = [line for line in structure_data if
                    'HOH' not in line]
    new_structure_data = protein_data + dowser_data
    return new_structure_data


def energy_minimize(dowser_data, structure_data):
    re_eval = open('re-eval.pdb', 'a')
    dowser_o = [line for line in dowser_data if ' OW ' in line]

    for i in tqdm(range(len(dowser_o))):
        current_water = dowser_o[i]
        with open('current_water.pdb', 'w') as cw:
            cw.write(current_water)
        water_pos_in_structure = [line for line, x in enumerate(structure_data)
                                  if x == current_water]
        water_pos_in_structure = water_pos_in_structure[0]
        current_structure = copy.deepcopy(structure_data)
        current_structure = water_hetatm_replacement(water_pos_in_structure,
                                                     current_structure)
        current_water_hetatm_OW = current_structure[water_pos_in_structure]
        current_water_hetatm_HW1 = current_structure[water_pos_in_structure+1]
        current_water_hetatm_HW2 = current_structure[water_pos_in_structure+2]
        print(f"water {i+1} is being changed to HETATM...\n", end='\r')
        print(current_water_hetatm_OW + current_water_hetatm_HW1
              + current_water_hetatm_HW2)
        with open('current_structure.pdb', 'w') as cs:
            cs.writelines(current_structure)

        run_reform('current_structure.pdb', 'current_structure_DOWSER.pdb')
        run_placeWat('current_structure_DOWSER.pdb',
                     'current_water.pdb', re_eval)

    refined_results = remove_disqualified_water('re-eval.pdb', cutoff)
    updated_structure_data = update_structure_data(refined_results,
                                                   structure_data)
    os.remove('current_water.pdb')
    os.remove('current_structure.pdb')
    os.remove('current_structure_DOWSER.pdb')
    re_eval.close()
    os.remove('re-eval.pdb')
    return refined_results, updated_structure_data


if __name__ == "__main__":
    try:
        args = parser.parse_args()
        dowser_o_input = args.dowser_input
        structure_input = args.structure_input
        refined_pdb = args.refined_pdb
        cutoff = args.cutoff
        if not any(vars(args).values()):
            raise Exception
    except Exception:
        print('Usage: python refine_results.py -d dowser_input' +
              ' -s structure_input -p refined_pdb -c cutoff')
    try:
        os.remove('re-eval.pdb')
        os.remove(refined_pdb)
    except OSError:
        print("Did't find previous results")
    with open(dowser_o_input, 'r') as dowser_o:
        # dowser_data = read_dowser_water(dowser_o)
        dowser_data = [line for line in dowser_o.readlines()
                       if 'ATOM' and ' HOH ' in line]
        num_of_water = int(len(dowser_data) / 3)

    with open(structure_input, 'r') as structure:
        structure_data = structure.readlines()
    print(f"refining the energies of {num_of_water} water molecules...")
    i = 1
    while True:
        print(f"round {i}...")
        old_water_count = int(len(dowser_data) / 3)
        refined_results, structure_data = energy_minimize(dowser_data,
                                                          structure_data)
        dowser_data = refined_results
        new_water_count = int(len(refined_results) / 3)
        if old_water_count == new_water_count:
            break
        i += 1

    with open(refined_pdb, 'w') as rp:
        rp.writelines(refined_results)
    pass
