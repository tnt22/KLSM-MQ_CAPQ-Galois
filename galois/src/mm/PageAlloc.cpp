/** Page Allocator Implementation -*- C++ -*-
 * @file
 * @section License
 *
 * Galois, a framework to exploit amorphous data-parallelism in irregular
 * programs.
 *
 * Copyright (C) 2011, The University of Texas at Austin. All rights reserved.
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
 * @section Description
 *
 * @author Andrew Lenharth <andrewl@lenharth.org>
 */

#include "Galois/Runtime/mm/Mem.h"
#include "Galois/Runtime/ll/gio.h"
#include "Galois/Runtime/ll/StaticInstance.h"
#include "Galois/Threads.h"

#include <sys/mman.h>
#include <map>
#include <vector>
#include <numeric>

#if defined(GALOIS_USE_NUMA) && !defined(GALOIS_FORCE_NO_NUMA)
#define USE_NUMA
#endif

#ifdef USE_NUMA
#include <numa.h>
#endif

// mmap flags
static const int _PROT = PROT_READ | PROT_WRITE;
static const int _MAP_BASE = MAP_ANONYMOUS | MAP_PRIVATE;
#ifdef MAP_POPULATE
static const int _MAP_POP  = MAP_POPULATE | _MAP_BASE;
#endif
#ifdef MAP_HUGETLB
static const int _MAP_HUGE_POP = _MAP_BASE | MAP_HUGETLB | _MAP_POP;
static const int _MAP_HUGE = _MAP_BASE | MAP_HUGETLB;
#endif

namespace {
struct FreeNode {
  FreeNode* next;
};
 
typedef Galois::Runtime::LL::PtrLock<FreeNode, true> HeadPtr;
typedef Galois::Runtime::LL::CacheLineStorage<HeadPtr> HeadPtrStorage;

// Tracks pages allocated
struct PAState {
  std::vector<int> counts;
  std::map<void*, HeadPtr*> ownerMap;
  PAState() { 
    counts.resize(Galois::Runtime::LL::getMaxThreads(), 0);
  }
};

static Galois::Runtime::LL::StaticInstance<PAState> PA;

#ifdef __linux__
#define DoAllocLock true
#else
#define DoAllocLock false
#endif
static Galois::Runtime::LL::SimpleLock<DoAllocLock> allocLock;
static Galois::Runtime::LL::SimpleLock<true> dataLock;
static __thread HeadPtr* head = 0;

void allocFromOS(int numPages) {

  HeadPtr*& h = head;
  if (!h) { //first allocation
     h = &((new HeadPtrStorage())->data);
  }

  allocLock.lock();
  dataLock.lock();
  h->lock();

  void* ptr = 0;

  for (int i = 0; i < numPages; i++) {
#ifdef MAP_HUGETLB
    //First try huge
    ptr = mmap(0, Galois::Runtime::MM::pageSize, _PROT, _MAP_HUGE_POP, -1, 0);
#endif

    //FIXME: improve failure case to ensure pageSize alignment
#ifdef MAP_POPULATE
    //Then try populate
    if (!ptr || ptr == MAP_FAILED)
      ptr = mmap(0, Galois::Runtime::MM::pageSize, _PROT, _MAP_POP, -1, 0);
#endif
    //Then try normal
    if (!ptr || ptr == MAP_FAILED) {
      ptr = mmap(0, Galois::Runtime::MM::pageSize, _PROT, _MAP_BASE, -1, 0);
    }

    if (!ptr || ptr == MAP_FAILED) {
      GALOIS_SYS_DIE("Out of Memory");
    }

    //protect the tracking structures
    PAState& p = *PA.get();
    p.ownerMap[ptr] = h;
    p.counts[Galois::Runtime::LL::getTID()] += 1;

    FreeNode* nh = reinterpret_cast<FreeNode*>(ptr);
    nh->next = h->getValue();
    h->setValue(nh);
  }

  h->unlock();
  dataLock.unlock();
  allocLock.unlock();
}

} // end anon namespace

void Galois::Runtime::MM::pageIn(void* buf, size_t len) {
  volatile char* ptr = reinterpret_cast<volatile char*>(buf);
  for (size_t i = 0; i < len; i += smallPageSize)
    ptr[i];
}

void* Galois::Runtime::MM::pageAlloc() {
  HeadPtr* phead = head;

  if (!phead || !phead->getValue()) {
    allocFromOS(16);
    phead = head;
  }
  if (phead) {
    phead->lock();
    FreeNode* h = phead->getValue();
    if (h) {
      phead->unlock_and_set(h->next);
      return h;
    }
    phead->unlock();
  }
  abort();
}

void Galois::Runtime::MM::pageFree(void* m) {
  dataLock.lock();
  HeadPtr* phead = PA.get()->ownerMap[m];
  dataLock.unlock();
  assert(phead);
  phead->lock();
  FreeNode* nh = reinterpret_cast<FreeNode*>(m);
  nh->next = phead->getValue();
  phead->unlock_and_set(nh);
}

void Galois::Runtime::MM::pagePreAlloc(int numPages) {
  allocFromOS(numPages);
}

int Galois::Runtime::MM::numPageAllocTotal() {
  PAState& p = *PA.get();
  return std::accumulate(p.counts.begin(), p.counts.end(), 0);
}

int Galois::Runtime::MM::numPageAllocForThread(unsigned tid) {
  return PA.get()->counts[tid];
}

void* Galois::Runtime::MM::largeAlloc(size_t len, bool preFault) {
  size_t size = (len + pageSize - 1) & (~(size_t)(pageSize - 1));
  void * ptr = 0;

  allocLock.lock();
#ifdef MAP_HUGETLB
  ptr = mmap(0, size, _PROT, preFault ? _MAP_HUGE_POP : _MAP_HUGE, -1, 0);
# ifndef MAP_POPULATE
  if (ptr != MAP_FAILED && ptr && preFault) {
    pageIn(ptr, size); // XXX should use hugepage stride
  }
# endif
#endif
#ifdef MAP_POPULATE
  if (preFault && (!ptr || ptr == MAP_FAILED))
    ptr = mmap(0, size, _PROT, _MAP_POP, -1, 0);
#endif
  if (!ptr || ptr == MAP_FAILED) {
    ptr = mmap(0, size, _PROT, _MAP_BASE, -1, 0);
    if (ptr != MAP_FAILED && ptr && preFault) {
      pageIn(ptr, size);
    }
  }
#ifdef USE_NUMA
  bool interleave = true;
  char host[256];
  if (gethostname(host, sizeof(host)) < 0) {
    perror("gethostname");
    exit(1);
  }
  if (!strcmp(host,"tapuz40")) {
    interleave = Galois::Runtime::LL::getMaxPackageForThread(Galois::getActiveThreads() - 1) > 0;
  } else if (!strcmp(host, "pollux")) {
    interleave = true;
  }
  if (interleave)
    numa_interleave_memory(ptr, size, numa_all_nodes_ptr);
#else
#if !defined(GALOIS_FORCE_NO_NUMA)
  GALOIS_DIE("NUMA is disabled");
#endif
#endif
  allocLock.unlock();

  if (!ptr || ptr == MAP_FAILED)
    GALOIS_SYS_DIE("Out of Memory");
  return ptr;
}

void Galois::Runtime::MM::largeFree(void* m, size_t len) {
  size_t size = (len + pageSize - 1) & (~(size_t)(pageSize - 1));
  allocLock.lock();
  munmap(m, size);
  allocLock.unlock();
}
