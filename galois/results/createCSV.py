# get's a folder's name that contains all the basic test results and combines them into csv


import os
import sys
import subprocess

import numpy as py
import scipy as sp


def main(argv):

    lst_files = [argv[0]+"/"+file_name for file_name in os.listdir(argv[0])]
    cmd = ["python", "csv.py"]+lst_files
    out = subprocess.check_output(cmd)
    fout = open("results.csv", 'w')
    fout.write(out)


if __name__ == "__main__":
    main(sys.argv[1:])
