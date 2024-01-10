#!/bin/bash

DOWSER_O_INPUT=$1

STRUCTURE_INPUT=$2

NUM_OF_WATER=$(awk 'END {print NR}' $DOWSER_O_INPUT)

echo "refining the energies of $NUM_OF_WATER water molecules..."

for ((i=1; i<2; i++))
    do
    CURRENT_WATER=$(sed -n ${i}p $DOWSER_O_INPUT)
    echo "$CURRENT_WATER" > current_water.pdb
    CURRENT_WATER_HETATM=$(echo "$CURRENT_WATER" | sed s/ATOM\ \ /HETATM/1 )
    echo "$CURRENT_WATER_HETATM"
    sed 's/$CURRENT_WATER/$CURRENT_WATER_HETATM/g' $STRUCTURE_INPUT > current_structure.pdb
done