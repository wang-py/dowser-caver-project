from random import gauss
import numpy as np
import matplotlib.pyplot as plt

# script to reproduce the energy figure in the AutoDock paper

gauss1_w = -0.0356
gauss2_w = -0.00516
repul_w = 0.84
hydrophobic_w = -0.0351
hydrogen_bond_w = -0.587

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
    coeff = np.polyfit([0.5, 1.5],[1, 0],1)
    if d < 0.5:
        d = 1
    elif d > 1.5:
        d = 0
    else:
        d = coeff[0] * d + coeff[1]

    return hydrophobic_w * d

def hydrogen_bond(d):
    coeff = np.polyfit([-0.7, 0],[1, 0],1)
    if d < -0.7:
        d = 1
    elif d > 0:
        d = 0
    else:
        d = coeff[0] * d + coeff[1]

    return hydrogen_bond_w * d

def steric(d):
    sum = gauss1_w * gauss1(d) + gauss2_w * gauss2(d) + repul_w * repulsion(d)
    return sum

if __name__ == "__main__":
    num_pts = 101
    d = np.linspace(-1,6,num_pts)
    steric_energy = []
    steric_hydrophobic_energy = []
    steric_hbond_energy = []
    for i in range(num_pts):
        d_i = d[i]
        steric_energy.append(steric(d_i))
        steric_hydrophobic_energy.append(steric(d_i)+hydrophobic(d_i))
        steric_hbond_energy.append(steric(d_i)+hydrogen_bond(d_i))
    plt.plot(d, steric_energy, 'o', label="steric")
    plt.plot(d, steric_hydrophobic_energy, 'o', label="steric+hydrophobic")
    plt.plot(d, steric_hbond_energy, 'o', label="steric+h-bond")
    plt.ylim([-0.25, 0.1])
    plt.legend()
    plt.show()
    pass