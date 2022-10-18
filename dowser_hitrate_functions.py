import numpy as np
import sys

def dowser_water_pruning(dowser_water_pdb):
    """
    function that extracts the positions of oxygen atoms and gets rid of 
    duplicates in the input file
    ---------------------------------------------------------------------------
    dowser_water_pdb: str
    filename of the dowser water pdb
    ---------------------------------------------------------------------------
    Returns:
    dowser_water_arr: ndarray
    n x 4 array that contains the position and energy of the dowser predicted waters
    """
    with open(dowser_water_pdb, 'r') as dwp:
        data = dwp.readlines()
        data = np.array([line[32:] for line in data if ' O ' in line]).astype(str)
        # get rid of duplicates
        data = np.unique(data)
    
    dowser_water_xyz = np.array([line[0:23].split() for line in data]).astype(float)
    dowser_water_energy = np.array([line[28:-1] for line in data]).astype(float)
    dowser_water_arr = []
    for i in range(len(dowser_water_energy)):
        dowser_water_arr.append(np.append(dowser_water_xyz[i], dowser_water_energy[i]))
    
    dowser_water_arr =  np.array(dowser_water_arr).astype(float)

    return dowser_water_arr

def read_exp_water(exp_water_pdb):
    """
    reads the pdb file of experimental water and returns array of positions
    ----------------------------------------------------------------------------
    exp_water_pdb: str
    filename of experimental water pdb
    ----------------------------------------------------------------------------
    Returns:
    exp_water_arr: ndarray
    N x 3 array of experimental water positions
    """

    with open(exp_water_pdb, 'r') as ewp:
        exp_data = [line.split()[6:9] for line in ewp.readlines() if 'HETATM' in line]
    exp_data = np.array(exp_data)
    exp_water_arr = exp_data.astype(float)
    
    return exp_water_arr

def hit_detection(exp_water, dowser_water):
    pass

def calculate_distance(exp_water, dowser_water):
    pass

