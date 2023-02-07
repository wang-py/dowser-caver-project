from dowser_database_functions import *

if __name__ == "__main__":
    atomdict_file = sys.argv[1]
    gromos_file = sys.argv[2]
    atomdict_data = read_atomdict(atomdict_file)
    read_gromos_ff(gromos_file)
    pass