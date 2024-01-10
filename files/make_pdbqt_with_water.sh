#!/bin/bash

ending="qt"
pdbqt=$1$ending

if [[ ! -s $pdbqt ]]; then
 ./pythonsh prepare_receptor4.py -r $1 -A bonds_hydrogens -U nphs
fi
if [[ ! -s $pdbqt ]]; then
 ./pythonsh prepare_receptor4.py -r $1 -A hydrogens -U nphs
fi
if [[ ! -s $pdbqt ]]; then
 ./pythonsh prepare_receptor4.py -r $1 -U nphs
fi

./pdbqt $pdbqt
cp new.pdbqt $pdbqt
rm new.pdbqt

