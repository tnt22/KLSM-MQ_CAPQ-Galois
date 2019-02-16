/** FIFO worklist -*- C++ -*-
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
#ifndef GALOIS_WORKLIST_GLOBPQ_H
#define GALOIS_WORKLIST_GLOBPQ_H

//#include "tbb/concurrent_priority_queue.h"
#include "WLCompileCheck.h"


namespace Galois {
namespace WorkList {

template<typename T = int, class PQT = LockFreeSkipList<DummyComparer<int>,int>, bool Concurrent = true>
struct GlobPQ : private boost::noncopyable {
  template<bool _concurrent>
  struct rethread { typedef GlobPQ<T, PQT, _concurrent> type; };

  template<typename _T>
  struct retype { typedef GlobPQ<_T, PQT, Concurrent> type; };

private:
  PQT pq; // tentative vertices

public:
  typedef T value_type;

  void push(const value_type& val) {
    pq.push(val);
  }

  template<typename Iter>
  unsigned int push(Iter b, Iter e) {
    int npush;
    for (npush = 0; b != e; npush++)
      push(*b++);
    return npush;
  }

  template<typename RangeTy>
  unsigned int push_initial(const RangeTy& range) {
    auto rp = range.local_pair();
    return push(rp.first, rp.second);
  }

  Galois::optional<value_type> pop() {
    Galois::optional<value_type> retval;
    value_type r;
    
    if (pq.try_pop(r))
      retval = r;
    return retval;
  }

  GlobPQ() { }
  GlobPQ(int qid) { }

};
GALOIS_WLCOMPILECHECK(GlobPQ)


} // end namespace WorkList
} // end namespace Galois

#endif
