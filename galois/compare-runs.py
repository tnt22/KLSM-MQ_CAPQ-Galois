#!/usr/bin/python
import sys, getopt, commands, os, datetime, signal, subprocess, shlex, time

def usage(val):
  sys.exit(val)

PUSH = 0
POP = 1

s0 = {}
s1 = {}

def avg(l):
  return sum(l)/len(l)

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

  fn = args[0]
  for ln in open(fn).readlines():
    if "STAT," in ln:
      k = ln.split(",")[2]
      l = s0.setdefault(k, [])
      l.append(float(ln.split(",")[4]))

  fn = args[1]
  for ln in open(fn).readlines():
    if "STAT," in ln:
      k = ln.split(",")[2]
      l = s1.setdefault(k, [])
      l.append(float(ln.split(",")[4]))

  for k in s0.keys():
    if k not in s1:
      continue
    if sum(s1[k]) == 0:
      continue
    if "Time" in k or "Work" in k:
      print "%-10s: %.3f [%.2f vs %.2f]" % (k, avg(s0[k])/avg(s1[k]), avg(s0[k])/avg(s0["Work"]), avg(s1[k])/avg(s1["Work"]))
    else:
      print "%-10s: %.3f" % (k, avg(s0[k])/avg(s1[k]))

if __name__ == "__main__":
  main(sys.argv[1:])
