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
    num_of_pts = tp.shape[0]
    final_output_name = output_folder + "/AllPredictions.pdbqt"
    # start docking water molecules along the caver tunnel
    for i in range(num_of_pts):
        prediction_filename = output_folder + "/point_%s.pdbqt"%str(i)
        prediction_temp = []
        vina_args = ("./vina", "--config", output_folder + "/box_%s.txt"%str(i), \
            "--out", prediction_filename)
        popen = subprocess.Popen(vina_args, stdout=subprocess.PIPE)
        popen.wait()
        output = popen.stdout.read()
        print(output.decode("utf-8"))
        # grepping the results
        with open(prediction_filename, 'r') as pf:
            for line in pf:
                if "ATOM" in line:
                    prediction_temp.append(line)

        # concatenate all predictions into one file
        with open(final_output_name, 'a+') as fo:
            fo.writelines(prediction_temp)

    pass