import sys
import os

if __name__ == "__main__":
    dowser_o_input=sys.argv[1]
    structure_input=sys.argv[2]
    with open(dowser_o_input) as dowser_o:
        data = dowser_o.readlines()
        num_of_water=len(data)

    print(f"refining the energies of {num_of_water} water molecules...")

    for i in range(num_of_water):
        current_water=data[i]
        print(current_water)
        #CURRENT_WATER_HETATM=$(echo "$CURRENT_WATER" | sed s/ATOM\ \ /HETATM/1 )
        #echo "$CURRENT_WATER_HETATM"
        #sed 's/$CURRENT_WATER/$CURRENT_WATER_HETATM/g' $STRUCTURE_INPUT > current_structure.pdb