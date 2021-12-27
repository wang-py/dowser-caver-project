import numpy as np
import pandas as pd
import sys
import os

#TODO: a function that scans the tunnel profiles csv and chooses a certain 
#tunnel. This function should output an array of coordinates and sizes.

#TODO: a function that reads the output of the above function, then outputs
#a series of config files for AutoDock Vina. Default box size is 6, but it
#could be customized to each sphere's size.

def csv_to_tunnel_points_arrays(tunnel_csv, tunnel_index):
    # skips information before coordinates
    skip_cols = 13
    # only takes coordinates of that tunnel
    tunnel_index_in_csv = (tunnel_index - 1) * 7 + 1
    rows_to_keep = [x for x in range(tunnel_index_in_csv, tunnel_index_in_csv + 3)]
    # scanning through the csv file
    tunnel_points = pd.read_csv(tunnel_csv, skiprows = lambda x: x not in rows_to_keep, \
        header=None)
    tunnel_points = tunnel_points.to_numpy()

    return tunnel_points[:, skip_cols:].astype(float)

def tunnel_points_to_box_configs(tunnel_points, output_folder):
    num_of_pts = tunnel_points.shape[1]
    receptor = "files/cavity_subunits.pdbqt"
    ligand = "files/water.pdbqt"
    energy_range = 100
    size = 6
    exhaustiveness = 20
    for i in range(num_of_pts):
        tp = tunnel_points[:,i]
        with open(output_folder + "/box_" + str(i) + ".txt", 'w') as file:
            file.write("receptor = " + str(receptor) + "\n")
            file.write("ligand = " + str(ligand) + "\n")
            file.write("\n")
            file.write("center_x = %.3f\n"%tp[0])
            file.write("center_y = %.3f\n"%tp[1])
            file.write("center_z = %.3f\n"%tp[2])
            file.write("\n")
            file.write("size_x = %.3f\n"%size)
            file.write("size_y = %.3f\n"%size)
            file.write("size_z = %.3f\n"%size)
            file.write("\n")
            file.write("energy_range = " + str(energy_range) + "\n")
            file.write("exhaustiveness = " + str(exhaustiveness) + "\n")
    pass

if __name__ == "__main__":
    tunnel_csv = sys.argv[1]
    tunnel_index = int(sys.argv[2])
    output_folder = "boxes"
    try:
        os.mkdir(output_folder)
    except OSError as error:
        print("folder " + output_folder + " was already created.\n")
    tp = csv_to_tunnel_points_arrays(tunnel_csv, tunnel_index)
    tunnel_points_to_box_configs(tp, output_folder)
    pass