from random import gauss
import numpy as np
import matplotlib.pyplot as plt

# script to reproduce the energy figure in the AutoDock paper

gauss1_w = -0.0356
gauss2_w = -0.00516
repul_w = 0.84
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
    sum = gauss1_w * gauss1(d) + gauss2_w * gauss2(d) + repul_w * repulsion(d)
    return sum

if __name__ == "__main__":
    num_pts = 41
    d = np.linspace(-1,6,num_pts)
    steric_energy = []
    for i in range(num_pts):
        d_i = d[i]
        steric_energy.append(steric(d_i))

    steric_energy = np.array(steric_energy)
    plt.plot(d, steric_energy, 'o', label="steric")
    plt.ylim([-0.25, 0.1])
    #plt.plot(d, steric_hydrophobic_energy, 'o', label="steric+hydrophobic")
    #plt.plot(d, steric_hbond_energy, 'o', label="steric+h-bond")
    plt.legend()
    plt.show()
    pass