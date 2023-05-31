from caver_to_vina import *
import matplotlib.pyplot as plt
import sys
# Script that reads the csv output of caver and plots tunnel profiles

def plotting(x, y, label):
    """
    plotting function
    ---------------------------------------------------------------------------
    x: ndarray
    x axis of the plot

    y: ndarray
    y axis of the plot

    label: str
    label of the data
    """
    #plt.figure()
    fig, ax = plt.subplots()
    ax.scatter(x, y, label=label)
    ax.hlines(1.4, np.min(x), np.max(x), linestyles='dashed', colors='k',label="radius of water")
    ax.set_xlabel("points")
    ax.set_ylabel("Radius [A]")
    ax.set_title("sphere size along caver points")
    ax.legend()

if __name__ == "__main__":
    # input tunnel pdb file
    tunnel_pdb = sys.argv[1]
    # figure title
    #figure_title = sys.argv[2]
    tp = pdb_to_tunnel_points_arrays(tunnel_pdb)
    radii = tp[:, 3]
    site_number = np.arange(tp.shape[0]) + 1
    plotting(site_number, radii, "E channel")
    plt.show()

    pass