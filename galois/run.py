#!/usr/bin/python
"""
Runs Galois experiments.

Usage: python run.py -g graph [-options] benchmarks...

Options:
  -m ..., --malloc=...      malloc library to use (default: libc)
  -t ..., --threads=...     # of threads
  -r ..., --runs=...        # of runs per benchmarks
  -v ..., --variant=...     variant to use (appended to benchmark name)
  -d ..., --directory=...   directory for outputs (default: results)
  -g ..., --graph=...       input graph
  -n, --numa                NUMA run (don't fake single board topology)
"""
import sys, getopt, commands, os, datetime, signal, subprocess, shlex, time, platform

def usage(val):
  print __doc__
  sys.exit(val)

class Alarm(Exception):
  pass

def alarm_handler(signum, frame):
  raise Alarm

def run(cmd, output, timeout=None, **params):
  signal.signal(signal.SIGALRM, alarm_handler)
  if timeout is not None:
    signal.alarm(timeout)
  output.write("==== BENCHMARK START ===\n")

  s = ""
  for param, value in params.iteritems():
    s += "%s=%s " % (param, value)
  if s != "":
    output.write("==== PARAMS %s ===\n" % (s))
  output.flush()

  start = time.time()
  try:
    proc = subprocess.Popen(shlex.split(cmd), shell=False, stdout=output, stderr=subprocess.STDOUT)
    r, err = proc.communicate()
    signal.alarm(0)
    print r
    if err is not None: print err
  except Alarm:
    print "Command '%s' timed out" % (cmd)
    commands.getoutput("kill -9 %d" % (proc.pid))
    proc.wait()
    output.write("STAT SINGLE Time (null) TIMEOUT-%d\n" % (timeout))
  elapsed = (time.time() - start)
  output.write("==== BENCHMARK TOOK %.3f\n" % (elapsed))
  output.close()

# ../../apps/pmst/Pmst.cpp:#include "Exp/PriorityScheduling/WorkListTL.h"

# Galois benchmarks
benchmarks = {
  "pagerank" : [ "pagerank/pagerank -t %d -wl %s -algo async_prt -graphTranspose $GRAPH.transpose -amp 100 -tolerance 0.001 $GRAPH", None ],
  "sssp" : [ "sssp/sssp -noverify -t %d -delta %d -wl %s -startNode NODE $GRAPH", 300 ],
  "matching" : [ "matching/bipartite-mcm -noverify -t %d -wl %s -file $GRAPH 1000000 100000000 10000 0", 120 ],
  "boruvkamerge" : [ "boruvka/boruvka-merge -noverify -t %d -wl %s $GRAPH", 120 ],
  "boruvka" : [ "boruvka/boruvka -noverify -t %d -wl %s $GRAPH", 120 ],
  "bfs" : [ "bfs/bfs -noverify -t %d -algo async -wl %s -startNode NODE $GRAPH", 300 ],
  "bfsbarrier" : [ "bfs/bfs -noverify -t %d -wl %s $GRAPH", 120 ],
  "betcet" : [ "betweennesscentrality/betweennesscentrality-inner -noverify -t %d -wl %s $GRAPH", 120 ],
  "betcetout" : [ "betweennesscentrality/betweennesscentrality-outer -noverify -t %d -wl %s $GRAPH", 120 ],
  "gmetis" : [ "gmetis/gmetis -noverify -t %d -wl %s $GRAPH 256", 180 ],
}

def main(argv):
  try:
    opts, args = getopt.getopt(argv, "m:t:r:v:d:g:D:n", [ "malloc=", "threads=", "runs=", "variant=", "directory=", "graph=", "delta=", "numa"])
  except getopt.GetoptError:
    usage(1)

  dir = "/specific/disk1/home/adamx/synch-proc/Galois-PQ/build/gcc5/apps/"
  env = ""

  if "gcc5" in dir:
    os.environ["LD_LIBRARY_PATH"] = "/usr/local/stow/gcc-5.3.0/lib/gcc-5.3.0/lib64" + ":" + os.environ["LD_LIBRARY_PATH"]
    os.environ["LIBRARY_PATH"] = "/usr/local/stow/gcc-5.3.0/lib/gcc-5.3.0/lib64" + ":" + os.environ["LIBRARY_PATH"]

  tests = args
  if len(tests) == 0:
    tests = benchmarks.keys()

  delta = 8
  for opt, arg in opts:
    if opt in ("-D", "--delta"):
      delta = int(arg)

  outdir = "results"
  for opt, arg in opts:
    if opt in ("-d", "--directory"):
      outdir = arg

  threads=20
  for opt, arg in opts:
    if opt in ("-t", "--threads"):
      threads = int(arg)

  runs = 5
  for opt, arg in opts:
    if opt in ("-r", "--runs"):
      runs = int(arg)

  wl = "obim"
  for opt, arg in opts:
    if opt in ("-v", "--variant"):
      wl = "%s" % (arg)

  numa = ""
  for opt, arg in opts:
    if opt in ("-n", "--numa"):
      env += "GALOIS_DONT_FAKE_TOPO= "
      numa = ".numa"

  malloc = "libc"
  for opt, arg in opts:
    if opt in ("-m", "--malloc"):
      malloc = arg

  graph = None
  inp = ""
  for opt, arg in opts:
    if opt in ("-g", "--graph"):
      graph = arg
      inp = "." + os.path.splitext(os.path.basename(arg))[0]

  if not graph:
    usage(1)

  lib = ""
  if malloc != "libc":
    env += "LD_PRELOAD=/specific/disk1/home/adamx/malloc/lib64/lib%s.so " % (malloc)

  if "klsm" in wl and lib == "":
    env += "LD_PRELOAD=/specific/disk1/home/adamx/malloc/lib64/libjemalloc.so "

  hostname = platform.uname()[1]
  kernel = platform.uname()[2]

  for test in tests:
    print test
    cmdline, timeout = benchmarks[test]
    if "-delta %d" in cmdline:
      cmdline = dir + cmdline % (threads, delta, wl)
    else:
      cmdline = dir + cmdline % (threads, wl)

    if "$GRAPH" in cmdline:
      cmdline = cmdline.replace("$GRAPH", graph)

    if "-startNode NODE" in cmdline:
      s = 0
      if "rmat-large.gr" in graph: s = 1
      if "twitter40.gr" in graph: s = 12
      cmdline = cmdline.replace("-startNode NODE", "-startNode %d" % (s))

    output = "%s/%s-%s%s%s.d%d_%d_%s.txt" % (outdir, test, wl, inp, numa, delta, threads, malloc)

    print "%s %s-%s %d threads" % (str(datetime.datetime.now()), test, wl, threads)
    cmd = "/usr/bin/env %s %s" % (env, cmdline)
    for i in range(0, runs):
      #commands.getoutput("./chkps.sh >>idiot.txt")
      if threads > 1:
        run(cmd, open(output, "a"), timeout, prog=test+"-"+wl+numa, hostname=hostname, kernel=kernel, malloc=malloc, threads=threads, delta=delta)
      else:
        run(cmd, open(output, "a"), timeout=600, prog=test+"-"+wl+numa, hostname=hostname, kernel=kernel, malloc=malloc, threads=threads, delta=delta)
      #if os.path.exists("core"):
      #  print "CORE DUMP!!"
      #  sys.exit(1)

if __name__ == "__main__":
  main(sys.argv[1:])
