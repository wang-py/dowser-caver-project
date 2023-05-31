# functions that modify the dowser database files atomdict.db and atomparms.db

import os
import sys

resname_gmx_to_dowser = {
    'GLUH': 'GLH ', 
    'ASPH': 'ASH ', 
    'LYSN': 'LYN ', 
    'HISD': 'HSN ', 
    'ARGN': 'ARN '}

class aminoacid:
    def __init__(self, resname, atom_list, charge_list):
        self.resname = resname
        self.atom_list = []
        selfcharge_list = []

def rename_neutral_residues(input_pdb, output_pdb):
    with open(input_pdb, 'r') as input:
        box = input.readline()
        data = [line for line in input.readlines() if 'ATOM' or 'HETATM' in line]
    with open(output_pdb, 'w') as output: 
        # keep box information
        output.write(box)
        # search every atom
        for line in data:
            for gmx in resname_gmx_to_dowser:
                dowser_resname = resname_gmx_to_dowser[gmx]
                # check if neutral residue names exist
                if gmx in line:
                    line = line.replace(gmx, dowser_resname)
            output.write(line)
    pass        

def read_atomdict(atomdict):
    with open(atomdict, 'r') as ad:
        atomdict_data = ad.readlines()
    atomdict_atom = [line for line in atomdict_data if 'ATOM' in line]
    return atomdict_atom

def read_gromos_ff(gromos_ff):
    with open(gromos_ff, 'r') as g_ff:
        g_ff_data = g_ff.readlines()
    
    g_ff_data_no_comments = [line.rstrip('\n') for line in g_ff_data if ';' not in line]
    pass
