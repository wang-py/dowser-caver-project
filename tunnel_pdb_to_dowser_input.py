from caver_to_vina import *

if __name__ == "__main__":
    tunnel_pdb = sys.argv[1]
    receptor = sys.argv[2]
    output_folder = "boxes"
    try:
        os.mkdir(output_folder)
    except OSError as error:
        print("folder " + output_folder + " was already created.\n")
    # convert pdb input into numpy arrays
    tp = pdb_to_tunnel_points_arrays(tunnel_pdb)
    # generate input files for AutoDock Vina
    tunnel_points_to_box_configs(tp, receptor, output_folder)

    pass