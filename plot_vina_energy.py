import matplotlib.pyplot as plt
import sys
import os
from glob import glob

if __name__ == "__main__":
    prediction_path = sys.argv[1]
    for pathname in glob(prediction_path + "/point_*.pdbqt"):
        print(pathname)
