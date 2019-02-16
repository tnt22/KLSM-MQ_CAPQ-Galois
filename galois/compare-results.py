#!/usr/bin/python
import sys, getopt, commands, os, datetime, signal, subprocess, shlex, time
from scipy.stats import gmean

def usage(val):
  sys.exit(val)

alld = []

def main(argv):
  try:
    opts, args = getopt.getopt(argv, "vr", ["validate", "rep"])
  except getopt.GetoptError:
    usage(1)

  validate = False
  for opt, arg in opts:
    if opt in ("-v", "--validate"):
      validate = True

  rep = False
  for opt, arg in opts:
    if opt in ("-r", "--rep"):
      rep = True

  dir = args[0]

  for fn in args[1:]:
    time0 = 0
    n0 = 0
    valid = True

    for ln in open(fn).readlines():
      if "==== BENCHMARK START ===" in ln:
        valid = True
      if "INVALID RUN?" in ln:
        valid = False
      if valid and "STAT," in ln and ",Time," in ln:
        time0 += float(ln.split(",")[4])
        n0 += 1

    if n0 == 0: continue

    time1 = 0
    n1 = 0

    if not os.path.isfile(dir+"/"+fn): continue

    for ln in open(dir+"/"+fn).readlines():
      if "==== BENCHMARK START ===" in ln:
        valid = True
      if "INVALID RUN?" in ln:
        valid = False
      if valid and "STAT," in ln and ",Time," in ln:
        time1 += float(ln.split(",")[4])
        n1 += 1

    if n1 == 0: continue

    time0 /= n0
    time1 /= n1
    d = time0/time1

    if (abs(1.0-d) > 0.05):
      if d <= 1:
        s = "BETTER"
      else:
        s = "WORSE"
      alld.append(d)
      print "%s now/%s %s -> %.3f" % (fn, dir, s, d)

  print "now/%s GEOMEAN %.3f" % (dir, gmean(alld))

if __name__ == "__main__":
  main(sys.argv[1:])
