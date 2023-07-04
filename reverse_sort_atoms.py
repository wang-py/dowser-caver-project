import numpy as np
import sys

def read_pdb(input_pdb):
    with open(input_pdb, 'r') as pdb_file:
        atom_info = [line for line in pdb_file.readlines() if 'ATOM' in line]
    
    return np.array(atom_info)

def reverse_atom_order(atom_list):
    rev_list = atom_list[::-1]
    return rev_list

def extract_residue(resname, atom_info):
    residue_atoms = [x for x in atom_info if resname == x.split()[5]]
    return residue_atoms

if __name__ == "__main__":
    atom_info = read_pdb(sys.argv[1])
    output_pdb = sys.argv[2]
    residues_to_reverse = ["62", "47"]
    for res in residues_to_reverse:
        res_atoms = extract_residue(res, atom_info)
        rev_res_atoms = reverse_atom_order(res_atoms)
        with open(output_pdb + "_" + res + ".pdb", 'w') as output_file:
            output_file.writelines(rev_res_atoms)
    
    pass
            