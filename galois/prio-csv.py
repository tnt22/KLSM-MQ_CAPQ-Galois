#!/usr/bin/python
import sys, getopt, commands, os, datetime, signal, subprocess, shlex, time

def usage(val):
  sys.exit(val)

PUSH = 0
POP = 1

def main(argv):
  try:
    opts, args = getopt.getopt(argv, "v", ["validate"])
  except getopt.GetoptError:
    usage(1)

  validate = False
  for opt, arg in opts:
    if opt in ("-v", "--validate"):
      validate = True

  if not validate:
    print "prog,malloc,threads,time,prio,nqueues,avg_qsize,npop,qpopfast,qpopfastcyc,qpoplocal,qpoplocalcyc,qpopremote,qpopremotecyc,qpopempty,qpopemptycyc"

  for fn in args:
    f = fn.split("_")
    prog = f[0].split(".")
    prog = "%s_%s_%s" % (prog[0], "-".join(prog[1:-1]), prog[-1])
    threads = int(f[1])
    malloc = os.path.splitext("-".join(f[2:]))[0]

    time = npop = 0
    qpopfast = {}
    qpopfastcyc = {}
    qpoplocal = {}
    qpoplocalcyc = {}
    qpopremote = {}
    qpopremotecyc = {}
    qpopempty = {}
    qpopemptycyc = {}
    timeout = False

    for ln in open(fn).readlines():
      if "STAT SINGLE Time (null) TIMEOUT" in ln:
        timeout = True
        continue
      if "STAT," in ln and ",Time," in ln:
        time = int(ln.split(",")[4])
      if "STAT," in ln and ",qPopEmpty," in ln:
        p = int(ln.split(",")[1])
        qpopempty[p] = float(ln.split(",")[4])
      if "STAT," in ln and ",qPopEmptyCyc," in ln:
        p = int(ln.split(",")[1])
        qpopemptycyc[p] = float(ln.split(",")[4])
      if "STAT," in ln and ",qPopFast," in ln:
        p = int(ln.split(",")[1])
        qpopfast[p] = float(ln.split(",")[4])
      if "STAT," in ln and ",qPopFastCyc," in ln:
        p = int(ln.split(",")[1])
        qpopfastcyc[p] = float(ln.split(",")[4])
      if "STAT," in ln and ",qPopLocal," in ln:
        p = int(ln.split(",")[1])
        qpoplocal[p] = float(ln.split(",")[4])
      if "STAT," in ln and ",qPopLocalCyc," in ln:
        p = int(ln.split(",")[1])
        qpoplocalcyc[p] = float(ln.split(",")[4])
      if "STAT," in ln and ",qPopRemote," in ln:
        p = int(ln.split(",")[1])
        qpopremote[p] = float(ln.split(",")[4])
      if "STAT," in ln and ",qPopRemoteCyc," in ln:
        p = int(ln.split(",")[1])
        qpopremotecyc[p] = float(ln.split(",")[4])
      if "==== BENCHMARK TOOK" in ln:
        nthreads = threads

        # have results?
        if time != 0:
          time /= 2.4*1000000000  # cycles to seconds

          if not validate:
            prios = qpopfast.keys()
            prios.sort()
            for i in xrange(1,len(prios)):
              this = prios[i]
              prev = prios[i-1]
              qpopfast[this] += qpopfast[prev]
              qpopfastcyc[this] += qpopfastcyc[prev]
              qpoplocal[this] += qpoplocal[prev]
              qpoplocalcyc[this] += qpoplocalcyc[prev]
              qpopremote[this] += qpopremote[prev]
              qpopremotecyc[this] += qpopremotecyc[prev]
              qpopempty[this] += qpopempty[prev]
              qpopemptycyc[this] += qpopemptycyc[prev]

            last = prios[-1]
            npop = qpopfast[last] + qpoplocal[last] + qpopremote[last]
            i = 1.0
            np = oldnp = oldn = 0
            for p in prios:
              np += 1
              n = qpopfast[p] + qpoplocal[p] + qpopremote[p]
              if (len(prios) < 1000) or (n >= (i/1000)*npop):
                frac = float(n)/npop
                print "%s,%s,%d,%s,%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (prog, malloc, threads, time, frac, np, float(n-oldn)/(np-oldnp), n, qpopfast[p], qpopfastcyc[p], qpoplocal[p], qpoplocalcyc[p], qpopremote[p], qpopremotecyc[p], qpopempty[p], qpopemptycyc[p])
                i = float(int(frac * 1000) + 1)
                oldnp = np
                oldn = n

        time = npop = 0
        qpopfast = {}
        qpopfastcyc = {}
        qpoplocal = {}
        qpoplocalcyc = {}
        qpopremote = {}
        qpopremotecyc = {}
        qpopempty = {}
        qpopemptycyc = {}
        timeout = False
        timeout = False

if __name__ == "__main__":
  main(sys.argv[1:])
