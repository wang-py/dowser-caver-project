# functions that modify the dowser database files atomdict.db and atomparms.db

import os
import sys

class aminoacid:
    def __init__(self, resname, atom_list, charge_list):
        self.resname = resname
        self.atom_list = []
        selfcharge_list = []

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