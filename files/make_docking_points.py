import numpy as np
import argparse

parser = argparse.ArgumentParser(prog="make_docking_points.py",
                                 description="script that generates" +
                                 " docking points based on protein structure")
parser.add_argument('-i', '--input_pdb', type=str)
parser.add_argument('-b', '--box_size', type=float, default=10.0)


def read_pdb(input_pdb):
    """
    reads a pdb file and returns numpy array of water data and protein data
    ----------------------------------------------------------------------------
    input_pdb: str
    path to pdb file
    ----------------------------------------------------------------------------
    Returns:
    atom_info: list
    python list that contains atom information from input pdb
    """
    # read in the pdb file
    pdb_file = open(input_pdb)
    atom_info = [line for line in pdb_file.readlines()
                 if line.startswith('ATOM  ')]

    return atom_info


def format_atom_info(atom_info):
    """
    reads a pdb file and returns numpy array of protein xyz
    ----------------------------------------------------------------------------
    atom_info: list
    python list that contains atom information from input pdb
    ----------------------------------------------------------------------------
    Returns:
    protein_data: ndarray: N x 3
    """
    protein_data = []
    for line in atom_info:
        one_data = np.array([])
        xyz = [float(x) for x in line[30:53].split()]

        one_data = np.append(one_data, xyz)
        protein_data.append(one_data)

    return np.array(protein_data)


def find_atoms_within_box(x, step_x, y, step_y, z, step_z, atoms):
    scan_x = np.logical_and(atoms[:, -3] < (x + step_x), atoms[:, -3] >= x)
    scan_y = np.logical_and(atoms[:, -2] < (y + step_y), atoms[:, -2] >= y)
    scan_z = np.logical_and(atoms[:, -1] < (z + step_z), atoms[:, -1] >= z)
    atoms_within_box = np.logical_and(scan_x, scan_y)
    atoms_within_box = np.logical_and(atoms_within_box, scan_z)

    return atoms_within_box


def get_input_partitions(atoms, partitions=2):
    """
    calculate partitions for protein atoms; the size is defined by box_size
    ----------------------------------------------------------------------------
    atoms: ndarray N x 3
    Array of protein atoms' information

    box_size: float angstoms
    Size of a single partition
    ----------------------------------------------------------------------------
    Returns:
    atoms_partitions: ndarray num_of_boxes x atoms_in_box x 7
    """
    box_min = np.min(atoms[:, -3:], axis=0)
    box_max = np.max(atoms[:, -3:], axis=0)
    # box_length = box_max - box_min
    # num_of_boxes = (box_length / box_size).astype(int)
    partitions_x, step_x = np.linspace(box_min[0], box_max[0], partitions,
                                       retstep=True)
    partitions_y, step_y = np.linspace(box_min[1], box_max[1], partitions,
                                       retstep=True)
    partitions_z, step_z = np.linspace(box_min[2], box_max[2], partitions,
                                       retstep=True)
    atoms_partitions = []
    for one_x in partitions_x:
        for one_y in partitions_y:
            for one_z in partitions_z:
                atoms_within_box = find_atoms_within_box(one_x, step_x,
                                                         one_y, step_y,
                                                         one_z, step_z,
                                                         atoms)
                if atoms_within_box.any():
                    atoms_partitions.append(([one_x, one_y, one_z],
                                             [one_x + step_x, one_y + step_y,
                                              one_z + step_z]))

    return atoms_partitions


def generate_one_docking_box(xyz, index: int, box_size: float = 10.0):
    box_prefix = "../boxes/box_"
    energy_range = 100
    exhaustiveness = 20
    num_modes = 20
    with open(box_prefix + f"{index}", 'w') as box:
        # file.write("receptor = " + str(receptor) + "\n")
        # file.write("ligand = " + str(ligand) + "\n")
        box.write("\n")
        box.write(f"center_x = {xyz[0]:.3f}\n")
        box.write(f"center_y = {xyz[1]:.3f}\n")
        box.write(f"center_z = {xyz[2]:.3f}\n")
        box.write("\n")
        box.write(f"size_x = {box_size}:.3f\n")
        box.write(f"size_y = {box_size}:.3f\n")
        box.write(f"size_z = {box_size}:.3f\n")
        box.write("\n")
        box.write("exhaustiveness = " + str(exhaustiveness) + "\n")
        box.write("energy_range = " + str(energy_range) + "\n")
        box.write("num_modes = " + str(num_modes) + "\n")
    pass


def generate_docking_box_in_one_partition(xyz_start, xyz_end, i: int,
                                          box_size: float = 10.0):
    boxes_x = np.arange(xyz_start[0], xyz_end[0], step=box_size)
    boxes_y = np.arange(xyz_start[1], xyz_end[1], step=box_size)
    boxes_z = np.arange(xyz_start[2], xyz_end[2], step=box_size)
    for one_x in boxes_x:
        for one_y in boxes_y:
            for one_z in boxes_z:
                generate_one_docking_box([one_x, one_y, one_z], index=i,
                                         box_size=box_size)

    pass


def generate_docking_boxes(input_partitions, box_size: float = 10.0):
    """
    generates docking box files based on input partitions
    ----------------------------------------------------------------------------
    input_partitions:
    list of tuples of bounds of input partitions

    box_size: float
    size of individual docking box
    ----------------------------------------------------------------------------
    Returns:
    txt files of docking boxes for AutoDock Vina
    """
    index = 1
    for p in input_partitions:
        generate_docking_box_in_one_partition(p[0], p[1], i=index,
                                              box_size=box_size)

    pass


if __name__ == "__main__":
    args = parser.parse_args()
    input_pdb = args.input_pdb
    box_size = args.box_size
    atom_info = read_pdb(input_pdb)
    atoms_xyz = format_atom_info(atom_info)
    input_partitions = get_input_partitions(atoms_xyz, partitions=4)

    pass
