/** Simple timer support -*- C++ -*-
 * @file
 * @section License
 *
 * Galois, a framework to exploit amorphous data-parallelism in irregular
 * programs.
 *
 * Copyright (C) 2013, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 *
 * @author Andrew Lenharth <andrewl@lenharth.org>
 */
#include "Galois/Timer.h"

#ifndef HAVE_CXX11_CHRONO
// This is linux/bsd specific
#include <sys/time.h>
#endif

using namespace Galois;

#ifndef HAVE_CXX11_CHRONO
Timer::Timer(bool on)
  :_start_hi(0), _start_low(0), _stop_hi(0), _stop_low(0)
{ if (on) start(); }

void Timer::start() {
  unsigned upper, lower;
  asm volatile ("rdtsc" : "=a"(lower), "=d"(upper));
  _start_hi = ((unsigned long)lower)|(((unsigned long)upper)<<32 );
}

void Timer::stop() {
  unsigned upper, lower;
  asm volatile ("rdtsc" : "=a"(lower), "=d"(upper));
  _stop_hi = ((unsigned long)lower)|(((unsigned long)upper)<<32 );
}

unsigned long Timer::sample() {
  stop();
  unsigned long t = get();
  _stop_hi = 0;
  return t;
}

unsigned long Timer::stopwatch() {
  stop();
  unsigned long t = get();
  _start_hi = _stop_hi;
  return t;
}
unsigned long Timer::get() const {
  return _stop_hi - _start_hi;
}

unsigned long Timer::get_usec() const {
  return _stop_hi - _start_hi;
}
#endif

TimeAccumulator::TimeAccumulator()
  :ltimer(), acc(0)
{}

void TimeAccumulator::start() {
  ltimer.start();
}

void TimeAccumulator::stop() {
  ltimer.stop();
  acc += ltimer.get_usec();
}

unsigned long TimeAccumulator::get() const {
  return acc;
}

TimeAccumulator& TimeAccumulator::operator+=(const TimeAccumulator& rhs) {
  acc += rhs.acc;
  return *this;
}

TimeAccumulator& TimeAccumulator::operator+=(const Timer& rhs) {
  acc += rhs.get_usec();
  return *this;
}
