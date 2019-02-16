#!/usr/bin/python
import sys, getopt, commands, os, datetime, signal, subprocess, shlex, time

def usage(val):
  sys.exit(val)

PUSH = 0
POP = 1

def main(argv):
  try:
    opts, args = getopt.getopt(argv, "vra", ["validate", "rep", "all"])
  except getopt.GetoptError:
    usage(1)

  validate = False
  printAll = False
  for opt, arg in opts:
    if opt in ("-v", "--validate"):
      validate = True
    if opt in ("-a", "--all"):
      printAll = True

  rep = False
  for opt, arg in opts:
    if opt in ("-r", "--rep"):
      rep = True

  if not validate:
    print "prog,malloc,threads,time,cycles,galois,push,pop,emptypop,npush,npop,nemptypop,good,bad,empty,work_cyc,ngood,nbad,nempty,work_n,qpopfast,qpopfastcyc,qpoplocal,qpoplocalcyc,qpopremote,qpopremotecyc,qpopempty,qpopemptycyc"

  for fn in args:
    f = fn.split("_")
    prog = f[0].split(".")
    prog = "%s_%s_%s" % (prog[0], "-".join(prog[1:-1]), prog[-1])
    threads = int(f[1])
    malloc = os.path.splitext("-".join(f[2:]))[0]

    time = cycles = galois = push = pop = emptypop = 0
    overall = good = bad = empty = work_cyc = 0
    noverall = ngood = nbad = nempty = work_n = 0
    npop = npush = nemptypop = 0
    qpopfast = qpopfastcyc = qpoplocal = qpoplocalcyc = qpopremote = qpopremotecyc = qpopempty = qpopemptycyc = 0
    conflicts = nconflicts = 0
    user = initpush = 0
    looptime = barrier = 0
    crit_barrier = {}
    crit_loop = {}
    crit_initpush = {}
    crit_work = {}
    perthread = {}
    timeout = False
    valid = True

    for ln in open(fn).readlines():
      if "STAT SINGLE Time (null) TIMEOUT" in ln:
        timeout = True
        continue
      if "possible co-runner" in ln and not printAll:
        valid = False
      if "STAT," in ln:
        loopname = ln.split(",")[1]
        loopvals = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",Time," in ln:
        time = float(ln.split(",")[4])
      if "STAT," in ln and ",LoopTime," in ln:
        looptime += float(ln.split(",")[4])
        crit_loop[loopname] = float(ln.split(",")[4])
      if "STAT," in ln and ",Work," in ln:
        cycles += float(ln.split(",")[4])
        crit_work[loopname] = loopvals[0]
        perthread["cycles"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",GaloisTime," in ln:
        galois += float(ln.split(",")[4])
        perthread["galois"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",PushTime," in ln:
        push += float(ln.split(",")[4])
        perthread["push"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",nPush," in ln:
        npush += float(ln.split(",")[4])
        perthread["npush"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",PopTime," in ln:
        pop += float(ln.split(",")[4])
        perthread["pop"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",nPop," in ln:
        npop += float(ln.split(",")[4])
        perthread["npop"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",EmptyPopTime," in ln:
        emptypop += float(ln.split(",")[4])
        perthread["emptypop"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and "nEmptyPop" in ln:
        nemptypop += float(ln.split(",")[4])
        perthread["nemptypop"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",EmptyWork," in ln:
        empty += float(ln.split(",")[4])
        perthread["empty"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",nEmpty," in ln:
        nempty += float(ln.split(",")[4])
        perthread["nempty"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",OverallWork," in ln:
        overall += float(ln.split(",")[4])
        perthread["overall"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",nOverall," in ln:
        noverall += float(ln.split(",")[4])
        perthread["noverall"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",BadWork," in ln:
        bad += float(ln.split(",")[4])
        perthread["bad"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",nBad," in ln:
        nbad += float(ln.split(",")[4])
        perthread["nbad"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",ConflictTime," in ln:
        conflicts += float(ln.split(",")[4])
        perthread["conflicts"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",Conflicts," in ln:
        nconflicts += float(ln.split(",")[4])
        perthread["nconflicts"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",UserTime," in ln:
        user += float(ln.split(",")[4])
        perthread["user"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",BarrierTime," in ln:
        barrier += float(ln.split(",")[4])
        crit_barrier[loopname] = loopvals[0]
        perthread["barrier"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",InitPushTime," in ln:
        initpush += float(ln.split(",")[4])
        crit_initpush[loopname] = loopvals[0]
        perthread["initpush"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",qPopEmpty," in ln:
        qpopempty += float(ln.split(",")[4])
        perthread["qpopempty"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",qPopEmptyCyc," in ln:
        qpopemptycyc += float(ln.split(",")[4])
        perthread["qpopemptycyc"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",qPopFast," in ln:
        qpopfast += float(ln.split(",")[4])
        perthread["qpopfast"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",qPopFastCyc," in ln:
        qpopfastcyc += float(ln.split(",")[4])
        perthread["qpopfastcyc"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",qPopLocal," in ln:
        qpoplocal += float(ln.split(",")[4])
        perthread["qpoplocal"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",qPopLocalCyc," in ln:
        qpoplocalcyc += float(ln.split(",")[4])
        perthread["qpoplocalcyc"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",qPopRemote," in ln:
        qpopremote += float(ln.split(",")[4])
        perthread["qpopremote"] = [float(x) for x in ln.split(",")[5:]]
      if "STAT," in ln and ",qPopRemoteCyc," in ln:
        qpopremotecyc += float(ln.split(",")[4])
        perthread["qpopremotecyc"] = [float(x) for x in ln.split(",")[5:]]
      if "==== BENCHMARK TOOK" in ln:
        # sanity check, since we changed accounting method
        if user != 0 and overall != 0:
          print "%s: have both user and overall" % (fn)
          sys.exit(1)
        if bad != 0 and conflicts != 0:
          print "%s: have both bad and conflicts" % (fn)
          sys.exit(1)

        nthreads = threads
        if time != 0 and rep:
          # find thread with median amount of work
          nthreads = 1
          m = [p for p in perthread["npop"] if p > 100]
          m.sort()
          m = m[len(m)/2]
          t = perthread["npop"].index(m)
          cycles = perthread["cycles"][t]
          galois = perthread["galois"][t]
          push = perthread["push"][t]
          npush = perthread["npush"][t]
          pop = perthread["pop"][t]
          npop = perthread["npop"][t]
          emptypop = perthread["emptypop"][t]
          nemptypop = perthread["nemptypop"][t]
          empty = perthread["empty"][t]
          nempty = perthread["nempty"][t]
          if overall != 0: overall = perthread["overall"][t]
          noverall = perthread["noverall"][t]
          bad = perthread["bad"][t]
          nbad = perthread["nbad"][t]
          conflicts = perthread["conflicts"][t]
          nconflicts = perthread["nconflicts"][t]
          if user != 0: user = perthread["user"][t]
          barrier = perthread["barrier"][t]
          initpush = perthread["initpush"][t]
          if qpopfast != 0:
            qpopempty = perthread["qpopempty"]
            qpopemptycyc = perthread["qpopemptycyc"]
            qpopfast = perthread["qpopfast"]
            qpopfastcyc = perthread["qpopfastcyc"]
            qpoplocal = perthread["qpoplocal"]
            qpoplocalcyc = perthread["qpoplocalcyc"]
            qpopremote = perthread["qpopremote"]
            qpopremotecyc = perthread["qpopremotecyc"]

        # have results?
        if time != 0:
          time = looptime

          if user == 0:
            good = overall - bad
            bad += conflicts
            work_cyc = good + bad + empty
          else:
            good = user - bad - empty
            bad += conflicts
            work_cyc = user + conflicts

          ngood = noverall - nbad
          nbad += nconflicts
          work_n = ngood + nbad + nempty

          # check that work is equal to the sum of its parts 
          if validate and noverall != 0 and cycles != 0:
            d = cycles/(work_cyc+galois+push+pop+emptypop)
            if (abs(1.0-d) > 0.1):
              print "work mismatch: %s,%.3f" % (fn,d)

          # handle initial pushes: add them to push time, and the rest
          # of the time, chalk up as galois time (this is time spent waiting
          # in the barriers due to load imbalance on initialization)
          if initpush != 0:
            # check that loop time equals to work cycles and initial pushes,
            # i.e., there are no "ghost cycles" spent elsewhere
            if validate and barrier:
              bad = None
              for l in crit_work:
                d = crit_loop[l] / (crit_work[l] + crit_initpush[l] + crit_barrier[l])
                if abs(1.0-d) > 0.1: bad = l
              if bad is not None:
                d = (looptime * threads) / (cycles + initpush + barrier)
                if abs(1.0-d) > 0.1:
                  print "missing looptime(%s): %s,%.3f" % (bad,fn,d)
            galois += (time * nthreads) - (work_cyc + push + pop + emptypop + galois + initpush)
            push += initpush

          # check that work and time count the same thing
          d = (work_cyc + galois + push + pop + emptypop) / (time * nthreads)
          if validate and noverall != 0 and (abs(1.0-d) > 0.1):
            print "work/time mismatch: %s,%.3f" % (fn,d)

          cycles = time * nthreads
          time /= 2.4*1000000000  # cycles to seconds
          if valid and not validate:
            print "%s,%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (prog, malloc, threads, time, cycles, galois, push, pop, emptypop, npush, npop, nemptypop, good, bad, empty, work_cyc, ngood, nbad, nempty, work_n, qpopfast, qpopfastcyc, qpoplocal, qpoplocalcyc, qpopremote, qpopremotecyc, qpopempty, qpopemptycyc)
        else:
            s = "NO-RESULTS"
            if timeout: s = "TIMEOUT"
            if validate: print "error (%s): %s" % (s, fn)

        time = cycles = galois = push = pop = emptypop = 0
        overall = good = bad = empty = work_cyc = 0
        noverall = ngood = nbad = nempty = work_n = 0
        npop = npush = nemptypop = 0
        qpopfast = qpopfastcyc = qpoplocal = qpoplocalcyc = qpopremote = qpopremotecyc = qpopempty = qpopemptycyc = 0
        conflicts = nconflicts = 0
        user = initpush = 0
        looptime = barrier = 0
        crit_barrier = {}
        crit_loop = {}
        crit_initpush = {}
        crit_work = {}
        timeout = False
        valid = True

if __name__ == "__main__":
  main(sys.argv[1:])
