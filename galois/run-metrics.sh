#!/bin/bash

#export LD_LIBRARY_PATH=/usr/lib/libblas:/usr/lib/atlas-base:${LD_LIBRARY_PATH}

results=results

rm corun.txt 

for g in scalefree/rmat16p-2e27.gr scalefree/rmat16p-2e24.gr # road/USA-road-d.USA.gr road/USA-road-t.USA.gr ljournal-2008.gr twitter40.gr random/r4-2e24.gr planar10M.bin scalefree/rmat-large.gr road/USA-road-d.USA.gr road/USA-road-t.USA.gr # clique4000.bin
do
  env OBIM_PRIO_STATS= ./run.py -m libc -t 1 -r 1 -g "inputs/${g}" -v obim -D 0 -d graph-metrics sssp bfs
done

