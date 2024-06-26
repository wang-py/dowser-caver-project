from caver_to_vina import *
import matplotlib.pyplot as plt
import sys
from pathlib import Path
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
    label_fontsize = 22
    fig, ax = plt.subplots(figsize=(10,7))
    ax.plot(x, y, label=label)
    ax.axhline(1.4, color='k', linestyle='--', label="radius of water")
    ax.tick_params(axis='x', labelsize=label_fontsize)
    ax.tick_params(axis='y', labelsize=label_fontsize)
    ax.set_xlabel("Caver channel points", fontsize=18)
    ax.set_ylabel("Radius [Å]", fontsize=18)
    #ax.set_title("sphere size along caver points", fontsize=label_fontsize)
    ax.legend(fontsize=label_fontsize)

if __name__ == "__main__":
    # input tunnel pdb file
    tunnel_pdb = sys.argv[1]
    input_filename = Path(tunnel_pdb).stem
    # figure title
    #figure_title = sys.argv[2]
    tp = pdb_to_tunnel_points_arrays(tunnel_pdb)
    radii = tp[:, 3]
    site_number = np.arange(tp.shape[0]) + 1
    plotting(site_number, radii, "E channel")
    plt.savefig(input_filename +"_paper.png", dpi=150)
    plt.show()

    pass