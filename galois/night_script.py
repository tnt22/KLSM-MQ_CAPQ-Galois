import os
import sys
import subprocess

import numpy as py
import scipy as sp


def main():
    """run the script python3 night_script.py"""

    folder = "testrun1"  # set the folder to save to, don't neeed to create one

    base = ["python", "run.py", "-g"]
    graphs = ["clique4000.gr", "ljournal-2008.gr", "planar10M.gr", "twitter40.gr",
              "random/r4-2e24.gr", "random/r4-2e26.gr", "road/USA-road-d.USA.gr",
              "road/USA-road-t.USA.gr", "scalefree/rmat16p-2e24.gr",
              "scalefree/rmat16p-2e27.gr"]
    klsm = ["klsm"+str(2**num) for num in range(8, 13)]
    fudim = ["obim"]+klsm
    kessler = ["capq"]
    leibo = ["multiqueue1", "multiqueue4"]
    data_structures = fudim  # choose one of the above
    applications = ["sssp"]
    delta = 4
    default = ["-m", "libc", "-D", str(delta), "-d", folder]
    threads = [str(num)
               for num in [1, 4, 8, 22, 44, 66, 88, 110, 132, 154, 176]]

    for graph in graphs:
        for data_structure in data_structures:
            for application in applications:
                for thread in threads:
                    files = os.listdir(folder)
                    graph_name = ""
                    if graph.find('/') != -1:
                        graph_name = graph[graph.find('/')+1:-3]
                    else:
                        graph_name = graph[:-3]

                    file_name = "{}-{}.{}.d{}_{}_{}.txt".format(
                        application, data_structure, graph_name, default[3], thread, default[1])
                    # print(file_name)
                    # print(str(files))
                    if(file_name not in files):
                        print(file_name)
                        # performs a test only if the file is not found in the folder
                        cmd = base+["/specific/disk1/home/mad/Galois-PQ/inputs/"+graph]+["-r", "1", "-t"]+[thread] + \
                            ["-v"]+[data_structure]+default+[application]
                        print("running: "+str(cmd))
                        out = subprocess.check_output(cmd)
                        files = os.listdir(folder)
                        paths = [os.path.join(folder, basename)
                                 for basename in files]
                        # returns the last created item from the directory
                        last = max(paths, key=os.path.getctime)
                        # print(last)
                        cmd = ["python", "csv.py", last]
                        out = subprocess.check_output(cmd)
                        # print(out.decode("utf-8"))
                        num_lines = len(out.decode("utf-8").split('\n'))
                        if(num_lines > 2):
                            # performs additonal 9 test if there wasn't a timeout in the first run
                            cmd = base+["/specific/disk1/home/mad/Galois-PQ/inputs/"+graph]+["-r", "9", "-t"]+[thread] + \
                                ["-v"]+[data_structure]+default+[application]
                            print("running: "+str(cmd))
                            out = subprocess.check_output(cmd)


if __name__ == "__main__":
    main()
