from dowser_database_functions import *

if __name__ == "__main__":
    gmx_pdb = sys.argv[1]
    dowser_pdb = sys.argv[2]
    renamed_data = rename_neutral_residues(gmx_pdb, dowser_pdb)