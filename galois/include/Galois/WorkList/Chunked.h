/** (d)Chunked(F|L)ifo worklist -*- C++ -*-
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

#ifndef GALOIS_WORKLIST_CHUNKED_H
#define GALOIS_WORKLIST_CHUNKED_H

#include "Galois/FixedSizeRing.h"
#include "Galois/Statistic.h"
#include "Galois/Timer.h"
#include "Galois/WorkList/WorkListHelpers.h"
#include "WLCompileCheck.h"

//#define PER_CHUNK_STATS

namespace Galois {
namespace WorkList {

//This overly complex specialization avoids a pointer indirection for non-distributed WL when accessing PerLevel
template<bool, template<typename> class PS, typename TQ>
struct squeue {
  PS<TQ> queues;
  TQ& get(int i) { return *queues.getRemote(i); }
  TQ& get() { return *queues.getLocal(); }
  int myEffectiveID() { return Runtime::LL::getTID(); }
  int size() { return Runtime::activeThreads; }
};

template<typename TQ>
struct squeue<true, Runtime::PerPackageStorage, TQ> {
  Runtime::PerPackageStorage<TQ> queues;
  TQ& get(int i) { return *queues.getRemote(Runtime::LL::getLeaderForPackage(i)); }
  TQ& get() { return *queues.getLocal(); }
  int myEffectiveID() { return Runtime::LL::getPackageForThread(Runtime::LL::getTID()); }
  int size() { return Runtime::LL::getMaxPackageForThread(Runtime::activeThreads-1) + 1; }
};

template<template<typename> class PS, typename TQ>
struct squeue<false, PS, TQ> {
  TQ queue;
  TQ& get(int i) { return queue; }
  TQ& get() { return queue; }
  int myEffectiveID() { return 0; }
  int size() { return 0; }
};

//! Common functionality to all chunked worklists
template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
struct ChunkedMaster : private boost::noncopyable {
  template<bool _concurrent>
  struct rethread { typedef ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, _concurrent> type; };

  template<typename _T>
  struct retype { typedef ChunkedMaster<_T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent> type; };

  template<int _chunk_size>
  struct with_chunk_size { typedef ChunkedMaster<T, QT, Distributed, DistStore, IsStack, _chunk_size, Concurrent> type; };

private:
  class Chunk : public FixedSizeRing<T, ChunkSize>, public QT<Chunk, Concurrent>::ListNode {};

  Runtime::MM::FixedSizeAllocator heap;

  struct p {
    Chunk* cur;
    Chunk* next;
    p(): cur(0), next(0) { }
  };

  typedef QT<Chunk, Concurrent> LevelItem;

  squeue<Concurrent, Runtime::PerThreadStorage, p> data;
  squeue<Distributed, DistStore, LevelItem> Q;

  Chunk* mkChunk() {
    return new (heap.allocate(sizeof(Chunk))) Chunk();
  }
  
  void delChunk(Chunk* C) {
    C->~Chunk();
    heap.deallocate(C);
  }

  void pushChunk(Chunk* C)  {
    LevelItem& I = Q.get();
    I.push(C);
  }

  Chunk* popChunkByID(unsigned int i)  {
    LevelItem& I = Q.get(i);
    return I.pop();
  }

  Chunk* popChunk(bool& local)  {
    int id = Q.myEffectiveID();
    Chunk* r = popChunkByID(id);
    if (r) {
      local = true;
      return r;
    }

    local = false;
    for (int i = id + 1; i < (int) Q.size(); ++i) {
      r = popChunkByID(i);
      if (r) 
	return r;
    }

    for (int i = 0; i < id; ++i) {
      r = popChunkByID(i);
      if (r)
	return r;
    }

    return 0;
  }

  template<typename... Args>
  T* emplacei(p& n, Args&&... args)  {
    T* retval = 0;
    if (n.next && (retval = n.next->emplace_back(std::forward<Args>(args)...)))
      return retval;
    if (n.next)
      pushChunk(n.next);
    Chunk* c = mkChunk();
    retval = c->emplace_back(std::forward<Args>(args)...);
    if (ChunkSize == 1) {
      pushChunk(c);
      assert(n.next == 0);
    } else {
      n.next = c;
    }
    assert(retval);
    return retval;
  }

public:
  typedef T value_type;

#ifdef PER_CHUNK_STATS
  Statistic *qPopFast, *qPopFastCyc, *qPopLocal, *qPopLocalCyc, *qPopRemote, *qPopRemoteCyc, *qEmpty, *qEmptyCyc;
#else
  static Statistic *qPopFast, *qPopFastCyc, *qPopLocal, *qPopLocalCyc, *qPopRemote, *qPopRemoteCyc, *qEmpty, *qEmptyCyc;
  static Runtime::LL::SimpleLock<true> statLock;
#endif

  void init_qstats(std::string id = "(NULL)") {
#ifndef PER_CHUNK_STATS
    if (qPopFast)
      return;

    statLock.lock();
    if (!qPopFast) {
#endif
      qPopFast = new Statistic("qPopFast", id);
      qPopFastCyc = new Statistic("qPopFastCyc", id);
      qPopLocal = new Statistic("qPopLocal", id);
      qPopLocalCyc = new Statistic("qPopLocalCyc", id);
      qPopRemote = new Statistic("qPopRemote", id);
      qPopRemoteCyc = new Statistic("qPopRemoteCyc", id);
      qEmpty = new Statistic("qPopEmpty", id);
      qEmptyCyc = new Statistic("qPopEmptyCyc", id);
#ifndef PER_CHUNK_STATS
    }
    statLock.unlock();
#endif
  }

  ChunkedMaster() : heap(sizeof(Chunk)) {
    init_qstats();
  }

  ChunkedMaster(int qid) : heap(sizeof(Chunk)) {
    if (getenv("OBIM_PRIO_STATS")) {
#ifndef PER_CHUNK_STATS
      GALOIS_DIE("need to define PER_CHUNK_STATS for OBIM_PRIO_STATS");
#endif
      init_qstats(std::to_string(qid));
    } else {
      init_qstats();
    }
  }

  ~ChunkedMaster() {
#ifndef PER_CHUNK_STATS
    if (!qPopFast)
      return;
    statLock.lock();
    if (qPopFast) {
#endif
      delete qPopFast;  qPopFast = 0;
      delete qPopFastCyc;  qPopFastCyc = 0;
      delete qPopLocal; qPopLocal = 0;
      delete qPopLocalCyc; qPopLocalCyc = 0;
      delete qPopRemote; qPopRemote = 0;
      delete qPopRemoteCyc; qPopRemoteCyc = 0;
      delete qEmpty; qEmpty = 0;
      delete qEmptyCyc; qEmptyCyc = 0;
#ifndef PER_CHUNK_STATS
    }
    statLock.unlock();
#endif
  }

  void flush() {
    p& n = data.get();
    if (n.next)
      pushChunk(n.next);
    n.next = 0;
  }
  
  /**
   * Construct an item on the worklist and return a pointer to its value.
   *
   * This pointer facilitates some internal runtime uses and is not designed
   * to be used by general clients. The address is generally not safe to use
   * in the presence of concurrent pops.
   */
  template<typename... Args>
  value_type* emplace(Args&&... args) {
    p& n = data.get();
    return emplacei(n, std::forward<Args>(args)...);
  }

  /**
   * Return pointer to next value to be returned by pop.
   *
   * For internal runtime use.
   */
  value_type* peek() {
    p& n = data.get();
    bool local;
    if (IsStack) {
      if (n.next && !n.next->empty())
	return &n.next->back();
      if (n.next)
	delChunk(n.next);
      n.next = popChunk(local);
      if (n.next && !n.next->empty())
	return &n.next->back();
      return NULL;
    } else {
      if (n.cur && !n.cur->empty())
	return &n.cur->front();
      if (n.cur)
	delChunk(n.cur);
      n.cur = popChunk(local);
      if (!n.cur) {
	n.cur = n.next;
	n.next = 0;
      }
      if (n.cur && !n.cur->empty())
	return &n.cur->front();
      return NULL;
    }
  }

  /**
   * Remove the value returned from peek() from the worklist. 
   *
   * For internal runtime use.
   */
  void pop_peeked() {
    p& n = data.get();
    if (IsStack) {
      n.next->pop_back();
      return;
    } else {
      n.cur->pop_front();
      return;
    }
  }

  void push(const value_type& val)  {
    p& n = data.get();
    emplacei(n, val);
  }

  template<typename Iter>
  unsigned int push(Iter b, Iter e) {
    p& n = data.get();
    int npush;
    for (npush = 0; b != e; npush++)
      emplacei(n, *b++);
    return npush;
  }

  template<typename RangeTy>
  unsigned int push_initial(const RangeTy& range) {
    auto rp = range.local_pair();
    return push(rp.first, rp.second);
  }

  Galois::optional<value_type> pop() {
    Timer tt(true);
    p& n = data.get();
    Galois::optional<value_type> retval;
    bool local;
    if (IsStack) {
      if (n.next && (retval = n.next->extract_back())) {
        *qPopFast += 1;
        *qPopFastCyc += tt.sample();
	return retval;
      }
      if (n.next)
	delChunk(n.next);
      n.next = popChunk(local);
      if (n.next) {
	retval = n.next->extract_back();
        if (local) {
          *qPopLocal += 1;
          *qPopLocalCyc += tt.sample();
        } else {
          *qPopRemote += 1;
          *qPopRemoteCyc += tt.sample();
        }
        return retval;
      }
      *qEmpty += 1;
      *qEmptyCyc += tt.sample();
      return Galois::optional<value_type>();
    } else {
      if (n.cur && (retval = n.cur->extract_front())) {
        *qPopFast += 1;
        *qPopFastCyc += tt.sample();
	return retval;
      }
      if (n.cur)
	delChunk(n.cur);
      n.cur = popChunk(local);
      if (!n.cur) {
	n.cur = n.next;
	n.next = 0;
      }
      if (n.cur) {
	retval = n.cur->extract_front();
        if (local) {
          *qPopLocal += 1;
          *qPopLocalCyc += tt.sample();
        } else {
          *qPopRemote += 1;
          *qPopRemoteCyc += tt.sample();
        }
        return retval;
      }
      *qEmpty += 1;
      *qEmptyCyc += tt.sample();
      return Galois::optional<value_type>();
    }
  }
};

#ifndef PER_CHUNK_STATS
template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent>::qPopFast;

template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent>::qPopFastCyc;

template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent>::qPopLocal;

template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent>::qPopLocalCyc;

template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent>::qPopRemote;

template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent>::qPopRemoteCyc;

template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent>::qEmpty;

template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
Runtime::LL::SimpleLock<true> ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent>::statLock;

template<typename T, template<typename, bool> class QT, bool Distributed, template<typename> class DistStore, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* ChunkedMaster<T, QT, Distributed, DistStore, IsStack, ChunkSize, Concurrent>::qEmptyCyc;
#endif

/**
 * Chunked FIFO. A global FIFO of chunks of some fixed size.
 *
 * @tparam ChunkSize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class ChunkedFIFO : public ChunkedMaster<T, ConExtLinkedQueue, false, Runtime::PerPackageStorage, false, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(ChunkedFIFO)

/**
 * Chunked LIFO. A global LIFO of chunks of some fixed size.
 *
 * @tparam ChunkSize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class ChunkedLIFO : public ChunkedMaster<T, ConExtLinkedStack, false, Runtime::PerPackageStorage, true, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(ChunkedLIFO)

/**
 * Distributed chunked FIFO. A more scalable version of {@link ChunkedFIFO}.
 *
 * @tparam ChunkSize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class dChunkedFIFO : public ChunkedMaster<T, ConExtLinkedQueue, true, Runtime::PerPackageStorage, false, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(dChunkedFIFO)

template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class dChunkedPTFIFO : public ChunkedMaster<T, ConExtLinkedQueue, true, Runtime::PerThreadStorage, false, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(dChunkedPTFIFO)

/**
 * Distributed chunked LIFO. A more scalable version of {@link ChunkedLIFO}.
 *
 * @tparam chunksize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class dChunkedLIFO : public ChunkedMaster<T, ConExtLinkedStack, true, Runtime::PerPackageStorage, true, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(dChunkedLIFO)

/**
 * Distributed chunked bag. A scalable and resource-efficient policy when you
 * are agnostic to the particular scheduling order.
 *
 * @tparam chunksize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class dChunkedBag : public ChunkedMaster<T, ConExtLinkedQueue, true, Runtime::PerPackageStorage, true, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(dChunkedBag)


} // end namespace WorkList
} // end namespace Galois

#endif
