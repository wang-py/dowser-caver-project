import sys
import os
import subprocess

def water_hetatm_replacement(position, structure_data):
    structure_data[position] = structure_data[position].replace('ATOM  ', 'HETATM')
    structure_data[position+1] = structure_data[position+1].replace('ATOM  ', 'HETATM')
    structure_data[position+2] = structure_data[position+2].replace('ATOM  ', 'HETATM')
    return structure_data

#./reform -pdbin lipolytica_JKAHN_charged_with_dowser_waters_HOH23_HETATM.pdb -pdbout lipolytica_JKAHN_charged_with_dowser_waters_HOH23_HETATM_DOWSER.pdb 
def run_reform(pdbin, pdbout):
    reform_args = ("./reform", "-pdbin", pdbin, "-pdbout", pdbout)
    popen = subprocess.Popen(reform_args, stdout=subprocess.PIPE)
    popen.wait()

#./placeWat lipolytica_JKAHN_charged_with_dowser_waters_HOH23_HETATM_DOWSER.pdb lipolytica_JKAHN_O23.pdb rotate > lipolytica_JKAHN_charged_with_dowser_waters_HOH23_HETATM_docking_1.pdb
def run_placeWat(dowser_pdb, current_water_pdb, output):
    placeWat_args = ('./placeWat', dowser_pdb, current_water_pdb, 'rotate')
    #print(placeWat_args)
    popen = subprocess.Popen(placeWat_args, stdout=output)
    popen.wait()

if __name__ == "__main__":
    dowser_o_input = sys.argv[1]
    structure_input = sys.argv[2]
    refined_pdb = sys.argv[3]
    with open(dowser_o_input, 'r') as dowser_o:
        dowser_data = dowser_o.readlines()
        num_of_water = len(dowser_data)
    
    with open(structure_input, 'r') as structure:
        structure_data = structure.readlines()

    refined = open(refined_pdb, 'a')

    print(f"refining the energies of {num_of_water} water molecules...")

    for i in range(1):
        current_water = dowser_data[i]
        print(current_water)
        with open('current_water.pdb', 'w') as cw:
            cw.write(current_water)
        water_pos_in_structure = [l for l,x in enumerate(structure_data) if x == current_water]
        water_pos_in_structure = water_pos_in_structure[0]
        current_structure = structure_data
        current_structure = water_hetatm_replacement(water_pos_in_structure, current_structure)
        current_water_hetatm_OW = current_structure[water_pos_in_structure]
        current_water_hetatm_HW1 = current_structure[water_pos_in_structure+1]
        current_water_hetatm_HW2 = current_structure[water_pos_in_structure+2]
        print(f"water {i+1} is being changed to HETATM...")
        print(current_water_hetatm_OW + current_water_hetatm_HW1 + current_water_hetatm_HW2)
        with open('current_structure.pdb', 'w') as cs:
            cs.writelines(current_structure)

        run_reform('current_structure.pdb', 'current_structure_DOWSER.pdb')
        run_placeWat('current_structure_DOWSER.pdb','current_water.pdb', refined)
