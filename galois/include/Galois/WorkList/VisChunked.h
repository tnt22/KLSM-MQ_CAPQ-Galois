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

#ifndef GALOIS_WORKLIST_VISCHUNKED_H
#define GALOIS_WORKLIST_VISCHUNKED_H

#include "Galois/FixedSizeRing.h"
#include "Galois/Statistic.h"
#include "Galois/Timer.h"
#include "Galois/Runtime/ll/PaddedLock.h"
#include "Galois/Runtime/ll/PtrLock.h"
#include "Galois/WorkList/WorkListHelpers.h"
#include "Galois/WorkList/Chunked.h"
#include "WLCompileCheck.h"

namespace Galois {
namespace WorkList {

//! Common functionality to all chunked worklists
template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
struct VisChunkedMaster : private boost::noncopyable {
  template<bool _concurrent>
  struct rethread { typedef VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, _concurrent> type; };

  template<typename _T>
  struct retype { typedef VisChunkedMaster<_T, QT, Distributed, IsStack, ChunkSize, Concurrent> type; };

  template<int _chunk_size>
  struct with_chunk_size { typedef VisChunkedMaster<T, QT, Distributed, IsStack, _chunk_size, Concurrent> type; };

private:
  class Chunk : public FixedSizeRing<T, ChunkSize>, public QT<Chunk, Concurrent>::ListNode {};

  Runtime::MM::FixedSizeAllocator heap;

  typedef Galois::Runtime::LL::PtrLock<Chunk, true> PtrLock;

  struct p {
    Chunk* cur;
    PtrLock next;
    p(): cur(0) { }
  };

  typedef QT<Chunk, Concurrent> LevelItem;

  squeue<Concurrent, Runtime::PerThreadStorage, p> data;
  squeue<Distributed, Runtime::PerPackageStorage, LevelItem> Q;

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
    PtrLock& L = n.next;
    L.lock();
    Chunk* next = L.getValue();
    if (next && (retval = next->emplace_back(std::forward<Args>(args)...))) {
      L.unlock();
      return retval;
    }
    if (next)
      pushChunk(next);
    next = mkChunk();
    retval = next->emplace_back(std::forward<Args>(args)...);
    L.unlock_and_set(next);
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

  VisChunkedMaster() : heap(sizeof(Chunk)) {
    init_qstats();
  }

  VisChunkedMaster(int qid) : heap(sizeof(Chunk)) {
    if (getenv("OBIM_PRIO_STATS")) {
#ifndef PER_CHUNK_STATS
      GALOIS_DIE("need to define PER_CHUNK_STATS for OBIM_PRIO_STATS");
#endif
      init_qstats(std::to_string(qid));
    } else {
      init_qstats();
    }
  }

  ~VisChunkedMaster() {
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
      PtrLock& L = n.next;
      L.lock();
      Chunk* next = L.getValue();
      if (next && (retval = next->extract_back())) {
        L.unlock();
        *qPopFast += 1;
        *qPopFastCyc += tt.sample();
	return retval;
      }
      if (next)
	    delChunk(next);
      next = popChunk(local);
      if (next) {
	    retval = next->extract_back();
        L.unlock_and_set(next);
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
        PtrLock& L = n.next;
        Chunk* next = L.getValue();
	    n.cur = next;
	    L.unlock_and_clear();
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
template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, Concurrent>::qPopFast;

template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, Concurrent>::qPopFastCyc;

template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, Concurrent>::qPopLocal;

template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, Concurrent>::qPopLocalCyc;

template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, Concurrent>::qPopRemote;

template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, Concurrent>::qPopRemoteCyc;

template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, Concurrent>::qEmpty;

template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
Runtime::LL::SimpleLock<true> VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, Concurrent>::statLock;

template<typename T, template<typename, bool> class QT, bool Distributed, bool IsStack, int ChunkSize, bool Concurrent>
Statistic* VisChunkedMaster<T, QT, Distributed, IsStack, ChunkSize, Concurrent>::qEmptyCyc;
#endif

/**
 * Chunked FIFO. A global FIFO of chunks of some fixed size.
 *
 * @tparam ChunkSize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class VisChunkedFIFO : public VisChunkedMaster<T, ConExtLinkedQueue, false, false, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(VisChunkedFIFO)

/**
 * Chunked LIFO. A global LIFO of chunks of some fixed size.
 *
 * @tparam ChunkSize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class VisChunkedLIFO : public VisChunkedMaster<T, ConExtLinkedStack, false, true, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(VisChunkedLIFO)

/**
 * Distributed chunked FIFO. A more scalable version of {@link ChunkedFIFO}.
 *
 * @tparam ChunkSize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class dVisChunkedFIFO : public VisChunkedMaster<T, ConExtLinkedQueue, true, false, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(dVisChunkedFIFO)

/**
 * Distributed chunked LIFO. A more scalable version of {@link ChunkedLIFO}.
 *
 * @tparam chunksize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class dVisChunkedLIFO : public VisChunkedMaster<T, ConExtLinkedStack, true, true, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(dVisChunkedLIFO)

/**
 * Distributed chunked bag. A scalable and resource-efficient policy when you
 * are agnostic to the particular scheduling order.
 *
 * @tparam chunksize chunk size
 */
template<int ChunkSize=64, typename T = int, bool Concurrent=true>
class dVisChunkedBag : public VisChunkedMaster<T, ConExtLinkedQueue, true, true, ChunkSize, Concurrent> {};
GALOIS_WLCOMPILECHECK(dVisChunkedBag)


} // end namespace WorkList
} // end namespace Galois

#endif
