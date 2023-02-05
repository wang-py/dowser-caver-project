# functions that modify the dowser database files atomdict.db and atomparms.db

import os
import sys

def read_atomdict(atomdict):
    with open(atomdict, 'r') as ad:
        atomdict_data = ad.readlines()
    atomdict_atom = [line for line in atomdict_data if 'ATOM' in line]
    return atomdict_atom

def read_gromos_ff(gromos_ff):
    pass