import numpy as np
import pandas as pd
import sys

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
    pass

if __name__ == "__main__":
    tunnel_csv = sys.argv[1]
    tunnel_index = int(sys.argv[2])
    output_folder = "boxes"
    tp = csv_to_tunnel_points_arrays(tunnel_csv, tunnel_index)
    tunnel_points_to_box_configs(tp, output_folder)
    pass