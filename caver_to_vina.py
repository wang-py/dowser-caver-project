import numpy as np
import pandas as pd
import sys

#TODO: a function that scans the tunnel profiles csv and chooses a certain 
#tunnel. This function should output an array of coordinates and sizes.

#TODO: a function that reads the output of the above function, then outputs
#a series of config files for AutoDock Vina. Default box size is 6, but it
#could be customized to each sphere's size.

def csv_to_tunnel_points_arrays(tunnel_csv):
    # skips information before coordinates
    skip_cols = 14
    # 2D array that stores all tunnel points' coordinates
    tunnel_points = []
    # temporary array for storing one point
    point_coords = []

    return tunnel_points

if __name__ == "__main__":
    tunnel_csv = sys.argv[1]
    pass