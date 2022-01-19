import numpy as np
import matplotlib.pyplot as plt

# script to reproduce the energy figure in the AutoDock paper

hydrophobic = -0.0351
hydrogen_bond = -0.587

def gauss1(d):
    pass

def gauss2(d):
    pass

def repulsion(d):
    pass

def scoring_func(d):
    pass

if __name__ == "__main__":
    d = np.linspace(0,2,21)
    energy = scoring_func(d)
    plt.plot(d, energy, 'o')
    plt.show()
    pass