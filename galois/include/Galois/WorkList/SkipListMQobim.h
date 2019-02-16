/** Scalable priority worklist -*- C++ -*-
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
#ifndef GALOIS_WORKLIST_SL_MQOBIM_H
#define GALOIS_WORKLIST_SL_MQBIM_H

#include "Galois/config.h"
#include "Galois/FlatMap.h"
#include "Galois/Timer.h"
#include "Galois/Statistic.h"
#include "Galois/Runtime/PerThreadStorage.h"
#include "Galois/WorkList/Fifo.h"
#include "Galois/WorkList/WorkListHelpers.h"

#include GALOIS_CXX11_STD_HEADER(type_traits)
#include <limits>

#include <iostream>

namespace Galois {
namespace WorkList {

/**
 * Approximate priority scheduling. Indexer is a default-constructable class
 * whose instances conform to <code>R r = indexer(item)</code> where R is
 * some type with a total order defined by <code>operator&lt;</code> and <code>operator==</code>
 * and item is an element from the Galois set iterator.
 *
 * An example:
 * \code
 * struct Item { int index; };
 *
 * struct Indexer {
 *   int operator()(Item i) const { return i.index; }
 * };
 *
 * typedef Galois::WorkList::OrderedByIntegerMetric<Indexer> WL;
 * Galois::for_each<WL>(items.begin(), items.end(), Fn);
 * \endcode
 *
 * @tparam Indexer Indexer class
 * @tparam Container Scheduler for each bucket
 * @tparam BlockPeriod Check for higher priority work every 2^BlockPeriod
 *                     iterations
 * @tparam BSP Use back-scan prevention
 */
template<class Indexer = DummyIndexer<int>, typename Container = FIFO<>,
  int BlockPeriod=0,
  typename T=int,
  typename Index=int,
  int C=4>
struct mqSkipListOrderedByIntegerMetric : private boost::noncopyable {
  template<bool _concurrent>
  struct rethread { typedef mqSkipListOrderedByIntegerMetric<Indexer, typename Container::template rethread<_concurrent>::type, BlockPeriod, T, Index> type; };

  template<typename _T>
  struct retype { typedef mqSkipListOrderedByIntegerMetric<Indexer, typename Container::template retype<_T>::type, BlockPeriod, _T, typename std::result_of<Indexer(_T)>::type> type; };

  template<unsigned _period>
  struct with_block_period { typedef mqSkipListOrderedByIntegerMetric<Indexer, Container, _period, T, Index> type; };

  template<typename _container>
  struct with_container { typedef mqSkipListOrderedByIntegerMetric<Indexer, _container, BlockPeriod, T, Index> type; };

  template<typename _indexer>
  struct with_indexer { typedef mqSkipListOrderedByIntegerMetric<_indexer, Container, BlockPeriod, T, Index> type; };

  typedef T value_type;

private:
  typedef typename Container::template rethread<false>::type CTy;
  typedef Galois::flat_map<Index, CTy*> LMapTy;
  //typedef std::map<Index, CTy*> LMapTy;

  struct Comparer: public std::binary_function<const int, const int, unsigned> {
    unsigned operator()(const int x, const int y) const {
      return x > y;
    }
  };

  struct perItem {
    LockFreeSkipListSet<Comparer, Index, CTy*> Q;
    Runtime::LL::PaddedLock<true> lock;
    Index curIndex;
    Index scanStart;
    CTy* current;
    unsigned int numPops;
    char pad[64];

    perItem() :
      curIndex(std::numeric_limits<Index>::min()), 
      scanStart(std::numeric_limits<Index>::min()),
      current(0), numPops(0) { }
  };

  typedef SkipListSetNode<Index,CTy*> sl_node_t;

  // NB: Place dynamically growing masterLog after fixed-size PerThreadStorage
  // members to give higher likelihood of reclaiming PerThreadStorage
  perItem* current;
  Galois::Timer clock;

  Runtime::MM::FixedSizeAllocator heap;
  Indexer indexer;
  int nQ;

  GALOIS_ATTRIBUTE_NOINLINE
  Galois::optional<T> globalPop() {
    Timer tt(true);
    unsigned tid = Galois::Runtime::LL::getTID();

    for (int i = 1; i < nQ; i++) {
      Galois::optional<T> rv;
      if (rv = slowPop(&current[(tid + i) % nQ])) {
        *GlobalPop += tt.sample();
        return rv;
      }
    }
    return Galois::optional<value_type>();
  }

  GALOIS_ATTRIBUTE_NOINLINE
  Galois::optional<T> randomPop() {
    Timer tt(true);

    *nRandomPop += 1;

    while (true) {
      int q0 = LockFreeSkipList<Comparer, Index>::rand_range(nQ) - 1;
      int q1 = LockFreeSkipList<Comparer, Index>::rand_range(nQ) - 1;
      perItem *p0 = &current[q0];
      perItem *p1 = &current[q1];
      Galois::optional<T> rv;
 
      *RandomPopIter += 1;

      if (p0->scanStart == std::numeric_limits<Index>::max() && 
          p1->scanStart == std::numeric_limits<Index>::max()) {
        break;
      } else if (p0->scanStart == std::numeric_limits<Index>::max()) {
        rv = slowPop(p1);
      } else if (p1->scanStart == std::numeric_limits<Index>::max()) {
        rv = slowPop(p0);
      } else if (p0->scanStart < p1->scanStart) {
        rv = slowPop(p0);
      } else {
        rv = slowPop(p1);
      }
      if (rv) {
        *RandomPop += tt.sample();
        return rv;
      }
    }
    return globalPop();
  }

  GALOIS_ATTRIBUTE_NOINLINE
  Galois::optional<T> slowPop(perItem* p) {

    Index msS = p->scanStart;
    if (msS == std::numeric_limits<Index>::max())
      return Galois::optional<value_type>();

    Timer tt(true);
    p->lock.lock();
    *popLock += tt.sample();

    sl_node_t *ii, *next;

    for (ii = p->Q.lower_bound(msS); ii->next[0]; ii = LockFreeSkipListSet<Comparer, Index, CTy*>::unset_mark(next)) {
      next = ii->next[0];
      if (LockFreeSkipListSet<Comparer, Index, CTy*>::is_marked(next))
        continue;

      Galois::optional<T> retval;
      if (ii->val && (retval = ii->val->pop())) {
        p->current = ii->val;
        p->curIndex = ii->key;
        p->scanStart = ii->key;
        tt.start();
        p->lock.unlock();
        *popLock += tt.sample();
        return retval;
      }
    }
    // nothing found, fail fast in future searches 
    p->scanStart = std::numeric_limits<Index>::max();
    tt.start();
    p->lock.unlock();
    *popLock += tt.sample();

    return Galois::optional<value_type>();
  }

  inline CTy* updateLocalOrCreate(perItem* p, Index i) {
    //Try local then try update then find again or else create and update the master log
    Timer tt(false);
    CTy* lC;
    tt.start();
    lC = p->Q.get(i);
    *readLocalCyc += tt.stopwatch();
    if (lC)
      return lC;
    lC = new (heap.allocate(sizeof(CTy))) CTy(i);
    *createCtyCyc += tt.stopwatch();

    tt.start();
    p->Q.push(i, lC);
    *updateLocalCyc += tt.stopwatch();

    return lC;
  }

  Statistic *updateLocalCyc, *readLocalCyc, *createCtyCyc, *popLock, *RandomPop, *GlobalPop, *RandomPopIter, *nRandomPop;

public:
  mqSkipListOrderedByIntegerMetric(const Indexer& x = Indexer()): heap(sizeof(CTy)), indexer(x), nQ(C*Galois::getActiveThreads()) {
    updateLocalCyc = new Statistic("ObimUpdateLocalCyc");
    readLocalCyc = new Statistic("ObimReadLocalCyc");
    createCtyCyc = new Statistic("ObimCreateCtyCyc");
    popLock = new Statistic("ObimPopLock");
    RandomPop = new Statistic("ObimRandomPopCyc");
    GlobalPop = new Statistic("ObimGlobalPopCyc");
    RandomPopIter = new Statistic("ObimRandomPopIter");
    nRandomPop = new Statistic("ObimNRandomPop");
    current = new perItem[nQ];
    clock.start();
  }

  ~mqSkipListOrderedByIntegerMetric() {
    sl_node_t *ii, *next;

    for (int i = 0; i < nQ; i++) {
      perItem *p = &current[i];
      for (ii = p->Q.head->next[0]; ii->next[0]; ii = LockFreeSkipListSet<Comparer, Index, CTy*>::unset_mark(next)) {
        next = ii->next[0];
        heap.deallocate(ii->val);
      }
    }
    delete updateLocalCyc;
    delete readLocalCyc;
    delete createCtyCyc;
    delete popLock;
    delete RandomPop;
    delete GlobalPop;
    delete RandomPopIter;
    delete nRandomPop;
  }

  void push(const value_type& val) {
    Index index = indexer(val);
    perItem* p;
    int q;

    do {
      q = LockFreeSkipList<Comparer, Index>::rand_range(nQ) - 1;
      p = &current[q];
    } while (!p->lock.try_lock());

    if (index == p->curIndex && p->current) {
      p->current->push(val);
      if (index < p->scanStart)
        p->scanStart = index;
      p->lock.unlock();
      return;
    }

    // Slow path
    CTy* lC = updateLocalOrCreate(p, index);
    if (index < p->scanStart)
      p->scanStart = index;
    // Opportunistically move to higher priority work
    if (index < p->curIndex) {
      p->curIndex = index;
      p->current = lC;
    }
    lC->push(val);
    p->lock.unlock();
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
    return randomPop();
  }
};
GALOIS_WLCOMPILECHECK(mqSkipListOrderedByIntegerMetric)

} // end namespace WorkList
} // end namespace Galois

#endif
