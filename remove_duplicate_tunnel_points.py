import numpy as np
import sys
from caver_to_vina import pdb_to_tunnel_points_arrays

def remove_duplicate_tp(tp_array):
    """
    removes duplicate tunnel points using numpy.unique() method
    """
    tp_unique = np.unique(tp_array, axis=0)

    return tp_unique

def format_atom_position(atom_pos):
    """
    takes an xyz array and formats it into pdb coordinates
    ----------------------------------------------------------------------------
    """
    x_pos = str("{0:.3f}".format(atom_pos[0])).rjust(8)
    y_pos = str("{0:.3f}".format(atom_pos[1])).rjust(8)
    z_pos = str("{0:.3f}".format(atom_pos[2])).rjust(8)
    atom_pos_pdb = x_pos + y_pos + z_pos

    return atom_pos_pdb

def write_edited_pdb(tp_array, output_pdb):
    """
    write edited pdb to file
    """
    with open(output_pdb, 'w') as output_file:
        for atom_i in range(tp_array.shape[0]):
            atom_i_pos = tp_array[atom_i, 0:3]
            atom_i_pos_pdb = format_atom_position(atom_i_pos)
            output_file.write("ATOM  {:>5}".format(atom_i+1) + "  H   FIL T   1    "
                              + atom_i_pos_pdb + "  0.00  " + "%.2f"%tp_array[atom_i, 3]
                              + "           H\n")
        output_file.write("END")

    pass

if __name__ == "__main__":
    tunnel_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    tp = pdb_to_tunnel_points_arrays(tunnel_pdb)
    tp_unique = remove_duplicate_tp(tp)
    write_edited_pdb(tp_unique, output_pdb)

    pass