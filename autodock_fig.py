from random import gauss
import numpy as np
import matplotlib.pyplot as plt

# script to reproduce the energy figure in the AutoDock paper

hydrophobic = -0.0351
hydrogen_bond = -0.587

def gauss1(d):
    return np.exp(-np.power(d/0.5, 2))

def gauss2(d):
    return np.exp(-np.power((d-3)/2, 2))

def repulsion(d):
    if d >= 0:
        d=0
    else:
        d=d*d
    
    return d

def hydrophobic(d):
    pass

def hydrogen_bond(d):
    pass

def steric(d):
    sum = gauss1(d) + gauss2(d) + repulsion(d)
    return sum

if __name__ == "__main__":
    d = np.linspace(0,2,21)
    steric_energy = steric(d)
    plt.plot(d, steric_energy, 'o', label="steric")
    #plt.plot(d, steric_hydrophobic_energy, 'o', label="steric+hydrophobic")
    #plt.plot(d, steric_hbond_energy, 'o', label="steric+h-bond")
    plt.legend()
    plt.show()
    pass