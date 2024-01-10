import sys
import os
import subprocess

if __name__ == "__main__":
    dowser_o_input=sys.argv[1]
    structure_input=sys.argv[2]
    with open(dowser_o_input, 'r') as dowser_o:
        dowser_data = dowser_o.readlines()
        num_of_water=len(dowser_data)
    
    with open(structure_input, 'r') as structure:
        structure_data = structure.readlines()

    print(f"refining the energies of {num_of_water} water molecules...")

    for i in range(1):
        current_water=dowser_data[i]
        print(current_water)
        with open('current_water.pdb', 'w') as cw:
            cw.write(current_water)
        #CURRENT_WATER_HETATM=$(echo "$CURRENT_WATER" | sed s/ATOM\ \ /HETATM/1 )
        #echo "$CURRENT_WATER_HETATM"
        #sed 's/$CURRENT_WATER/$CURRENT_WATER_HETATM/g' $STRUCTURE_INPUT > current_structure.pdb