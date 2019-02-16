/** Worklist building blocks -*- C++ -*-
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
#ifndef GALOIS_RUNTIME_WORKLISTHELPERS_H
#define GALOIS_RUNTIME_WORKLISTHELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <climits>
#include "WLCompileCheck.h"

#include "Galois/Runtime/Termination.h"
#include "Galois/Runtime/ll/PtrLock.h"

#include "k_lsm/k_lsm.h"

#include <boost/iterator/iterator_facade.hpp>
#include <boost/heap/d_ary_heap.hpp>

#define MEM_BARRIER     asm volatile("":::"memory")
#define ATOMIC_CAS_MB(p, o, n)  __sync_bool_compare_and_swap(p, o, n)
#define ATOMIC_FETCH_AND_INC_FULL(p) __sync_fetch_and_add(p, 1)

namespace Galois {
namespace WorkList {

template<typename T>
class ConExtListNode {
  T* next;
public:
  ConExtListNode() :next(0) {}
  T*& getNext() { return next; }
  T*const& getNext() const { return next; }
};

template<typename T>
class ConExtIterator: public boost::iterator_facade<
                      ConExtIterator<T>, T, boost::forward_traversal_tag> {
  friend class boost::iterator_core_access;
  T* at;

  template<typename OtherTy>
  bool equal(const ConExtIterator<OtherTy>& o) const { return at == o.at; }

  T& dereference() const { return *at; }
  void increment() { at = at->getNext(); }

public:
  ConExtIterator(): at(0) { }
  
  template<typename OtherTy>
  ConExtIterator(const ConExtIterator<OtherTy>& o): at(o.at) { }
  
  explicit ConExtIterator(T* x): at(x) { }
};

template<typename T, bool concurrent>
class ConExtLinkedStack {
  Runtime::LL::PtrLock<T, concurrent> head;
  
public:
  typedef ConExtListNode<T> ListNode;

  bool empty() const {
    return !head.getValue();
  }

  void push(T* C) {
    T* oldhead(0);
    do {
      oldhead = head.getValue();
      C->getNext() = oldhead;
    } while (!head.CAS(oldhead, C));
  }

  T* pop() {
    //lock free Fast path (empty)
    if (empty()) return 0;
    
    //Disable CAS
    head.lock();
    T* C = head.getValue();
    if (!C) {
      head.unlock();
      return 0;
    }
    head.unlock_and_set(C->getNext());
    C->getNext() = 0;
    return C;
  }

  //! iterators not safe with concurrent modifications
  typedef T value_type;
  typedef T& reference;
  typedef ConExtIterator<T> iterator;
  typedef ConExtIterator<const T> const_iterator;

  iterator begin() { return iterator(head.getValue()); }
  iterator end() { return iterator(); }

  const_iterator begin() const { return const_iterator(head.getValue()); }
  const_iterator end() const { return const_iterator(); }
};

template<typename T, bool concurrent>
class ConExtLinkedQueue {
  Runtime::LL::PtrLock<T,concurrent> head;
  T* tail;
  
public:
  typedef ConExtListNode<T> ListNode;
  
  ConExtLinkedQueue() :tail(0) { }

  bool empty() const {
    return !tail;
  }

  void push(T* C) {
    head.lock();
    //std::cerr << "in(" << C << ") ";
    C->getNext() = 0;
    if (tail) {
      tail->getNext() = C;
      tail = C;
      head.unlock();
    } else {
      assert(!head.getValue());
      tail = C;
      head.unlock_and_set(C);
    }
  }

  T* pop() {
    //lock free Fast path empty case
    if (empty()) return 0;

    head.lock();
    T* C = head.getValue();
    if (!C) {
      head.unlock();
      return 0;
    }
    if (tail == C) {
      tail = 0;
      assert(!C->getNext());
      head.unlock_and_clear();
    } else {
      head.unlock_and_set(C->getNext());
      C->getNext() = 0;
    }
    return C;
  }

  //! iterators not safe with concurrent modifications
  typedef T value_type;
  typedef T& reference;
  typedef ConExtIterator<T> iterator;
  typedef ConExtIterator<const T> const_iterator;

  iterator begin() { return iterator(head.getValue()); }
  iterator end() { return iterator(); }

  const_iterator begin() const { return const_iterator(head.getValue()); }
  const_iterator end() const { return const_iterator(); }
};

template<typename T>
struct DummyIndexer: public std::unary_function<const T&,unsigned> {
  unsigned operator()(const T& x) { return 0; }
};

template<typename T>
struct DummyComparer: public std::binary_function<const T&,const T&,unsigned> {
  unsigned operator()(const T& x, const T&y) { return x > y; }
};

// check cache alignment
template<typename K>
struct SkipListNode {
public:
  // dummy field which is used by the heap when the node is freed.
  // (without it, freeing a node would corrupt a field, possibly affecting
  // a concurrent traversal.)
  SkipListNode<K>* dummy;

  K key;
  int toplevel;
  SkipListNode* volatile next[0];

  void init(int _t, SkipListNode *_next) {
    toplevel = _t;
    for (int i = 0; i < _t; i++)
      next[i] = _next;
  }

  //T*& getNext() { return next; }
  //T*const& getNext() const { return next; }
};

#define SKIPLIST_LEVELS	24

template<class Comparer, typename K>
class LockFreeSkipList {

protected:
  typedef SkipListNode<K> sl_node_t;

  static Runtime::MM::ListNodeHeap heap[3];
  Runtime::TerminationDetection& term;
  Comparer compare;

  sl_node_t* head;
  uint8_t levelmax;

  static inline bool is_marked(sl_node_t* i)
  {
    return ((uintptr_t)i & (uintptr_t)0x01) != 0;
  }

  static inline bool is_dead(sl_node_t* i)
  {
    return ((uintptr_t)i & (uintptr_t)0x02) == 2;
  }

  static inline sl_node_t* unset_mark(sl_node_t* i)
  {
    return (sl_node_t *)((uintptr_t)i & ~(uintptr_t)0x03);
  }

  static inline sl_node_t* set_mark(sl_node_t* i)
  {
    return (sl_node_t *)((uintptr_t)i | (uintptr_t)0x01);
  }

  static inline sl_node_t * set_dead(sl_node_t * i)
  {
    return (sl_node_t *)((uintptr_t)i | (uintptr_t)0x02);
  }

  //Marsaglia's xorshf generator
  static inline unsigned long xorshf96(unsigned long* x, unsigned long* y, unsigned long* z)  //period 2^96-1
  {
    unsigned long t;
    (*x) ^= (*x) << 16;
    (*x) ^= (*x) >> 5;
    (*x) ^= (*x) << 1;

    t = *x;
    (*x) = *y;
    (*y) = *z;
    (*z) = t ^ (*x) ^ (*y);

    return *z;
  }

  static __thread unsigned long seeds[3];
  static __thread bool seeds_init;

public:
  static inline long rand_range(long r)
  {
    if (!seeds_init) {
      int fd = open("/dev/urandom", O_RDONLY);
      if (read(fd, seeds, 3 * sizeof(unsigned long)) < 0) {
        perror("read");
        exit(1);
      }
      close(fd); 
      seeds_init = true;
    }
    long v = xorshf96(seeds, seeds + 1, seeds + 2) % r;
    v++;
    return v;
  }

protected:
  void mark_node_ptrs(sl_node_t *n)
  {
    sl_node_t *n_next;
    int i;
    
    for (i=n->toplevel-1; i>=0; i--)
    {
      do
      {
        n_next = n->next[i];
        if (is_marked(n_next))
        {
          break;
        }
      } while (!ATOMIC_CAS_MB(&n->next[i], n_next, set_mark(n_next)));
    }
  }

  static __thread int spray_seed;
  static __thread bool spray_seed_init;

  static int _MarsagliaXOR(void) {

    if (!spray_seed_init) {
      int fd = open("/dev/urandom", O_RDONLY);
      if (read(fd, &spray_seed, sizeof(int)) < 0) {
        perror("read");
        exit(1);
      }
      close(fd);
      spray_seed_init = true;
    }

    const int a =      123456789;
    const int m =     2147483647;
    const int q =      521288629;  /* m div a */
    const int r =       88675123;  /* m mod a */
    int hi   = spray_seed / q;
    int lo   = spray_seed % q;
    int test = a * lo - r * hi;
    if (test > 0)
      spray_seed = test;
    else
      spray_seed = test + m;

    return spray_seed;
  }

  static int floor_log_2(unsigned int n) {
    int pos = 0;
    if (n >= 1<<16) { n >>= 16; pos += 16; }
    if (n >= 1<< 8) { n >>=  8; pos +=  8; }
    if (n >= 1<< 4) { n >>=  4; pos +=  4; }
    if (n >= 1<< 2) { n >>=  2; pos +=  2; }
    if (n >= 1<< 1) {           pos +=  1; }
    return ((n == 0) ? (-1) : pos);
  }

  inline sl_node_t *sl_new_node(int levelmax, sl_node_t *next) {
    int e = term.getEpoch() % 3;
    sl_node_t *node = reinterpret_cast<sl_node_t *>(heap[e].allocate(sizeof(sl_node_t) + levelmax*sizeof(sl_node_t*), levelmax-1));
    node->init(levelmax, next);
    return node;
  }

  inline sl_node_t *sl_new_node_key(K key, int levelmax) {
    sl_node_t *node = sl_new_node(levelmax, 0);
    node->key = key;
    return node;
  }

  inline void sl_delete_node(sl_node_t *n) {
    int e = (term.getEpoch() + 2) % 3;
    heap[e].deallocate(n, n->toplevel-1);
  }

public:

  LockFreeSkipList() : term(Runtime::getSystemTermination()), levelmax(SKIPLIST_LEVELS) {
    sl_node_t *min, *max;

    max = sl_new_node(levelmax, NULL);
    min = sl_new_node(levelmax, max);

    head = min;
  }

  bool empty() const
  {
    return head->next[0]->next[0] == 0;
  }

  void fraser_search(K key, sl_node_t **left_list, sl_node_t **right_list, sl_node_t *dead)
  {
    sl_node_t *left, *left_next, *right, *right_next;
    int i;

  retry:
    left = head;
    for (i = (dead ? dead->toplevel : levelmax) - 1; i >= 0; i--)
    {
      sl_node_t *first = NULL;

      left_next = left->next[i];
      if (is_marked(left_next))
        goto retry;
      /* Find unmarked node pair at this level */
      for (right = left_next; ; right = right_next)
      {
        /* Skip a sequence of marked nodes */
        while(1)
        {
          right_next = right->next[i];
          if (!is_marked(right_next))
            break;
          right = unset_mark(right_next);
        }
        /* Ensure left and right nodes are adjacent */
        if (left_next != right) {
          if (!ATOMIC_CAS_MB(&left->next[i], left_next, right))
            goto retry;
          for (sl_node_t *t = left_next; t != right; t = unset_mark(t->next[i]))
            t->next[i] = set_dead(t->next[i]);
        }
        /* When deleting, we have to keep going until we find our target node, or until
           we observe it has been deleted (right->key > key).  Once this happens, however,
           we need to descend to a node whose key is smaller than our target's, otherwise
           we might miss our target in the level below.  (Consider N1 and N2 with the same
           key, where on level 1, N1 -> N2 but on level 0, N2 -> N1; if when looking for N2
           we descend at N1, we'll miss N2 at level 0.) */
        if (!dead) {
          if (!right_next || !compare(key, right->key))
            break;
        } else {
          if (!first && !compare(key, right->key))
            first = left;
          if (!right_next || is_dead(dead->next[i]) || compare(right->key, key)) {
            if (first) left = first;
            break;
          }
        }
        left = right;
        left_next = right_next;
      }
      if (left_list != NULL)
        left_list[i] = left;
      if (right_list != NULL)
        right_list[i] = right;
    }
  }

  int get_rand_level()
  {
    int i, level = 1;
    for (i = 0; i < levelmax - 1; i++)
      {
        if ((rand_range(100)-1) < 50)
          level++;
        else
          break;
      }
    /* 1 <= level <= *levelmax */
    return level;
  }

  bool push(const K& key)
  {
    sl_node_t *newn, *new_next, *pred, *succ, *succs[levelmax], *preds[levelmax];
    int i, result = 0;

    newn = sl_new_node_key(key, get_rand_level());

  retry:
    fraser_search(key, preds, succs, NULL);
    if (succs[0]->key == key)
    {                             /* Value already in list */
      result = 0;
      sl_delete_node(newn);
      goto end;
    }

    for (i = 0; i < newn->toplevel; i++)
    {
      newn->next[i] = succs[i];
    }

    /* Node is visible once inserted at lowest level */
    if (!ATOMIC_CAS_MB(&preds[0]->next[0], succs[0], newn))
    {
      goto retry;
    }

    for (i = 1; i < newn->toplevel; i++)
    {
      while (1)
      {
        pred = preds[i];
        succ = succs[i];
        new_next = newn->next[i];
        /* Give up if pointer is marked */
        if (is_marked(new_next))
          goto success;
        /* Update the forward pointer if it is stale, which can happen
           if we called search again to update preds and succs. */
        if (new_next != succ && !ATOMIC_CAS_MB(&newn->next[i], new_next, succ))
          goto success;
        /* We retry the search if the CAS fails */
        if (ATOMIC_CAS_MB(&pred->next[i], succ, newn)) {
          if (is_marked(newn->next[i])) {
            fraser_search(key, NULL, NULL, newn);
            goto success;
          }
          break;
        }

        fraser_search(key, preds, succs, NULL);
      }
    }

  success:
    result = 1;

  end:
    return result;
  }

  sl_node_t* peek_pop(void) const
  {
    sl_node_t *first, *next;

    first = head;

    do {
      first = unset_mark(first->next[0]);
      next = first->next[0];
    } while(next && is_marked(next));

    return next ? first : 0;
  }

  K& get_min(void) const
  {
    sl_node_t *n = peek_pop();

    return n ? n->key : head->key;
  }

  bool complete_pop(sl_node_t *first, K& key)
  {
    sl_node_t *next = first->next[0];

    if (is_marked(next) ||
        !ATOMIC_CAS_MB(&first->next[0], next, set_mark(next)))
      return false;

    key = (first->key);
    mark_node_ptrs(first);

    fraser_search(key, NULL, NULL, first);
    sl_delete_node(first);

    return true;
  }

  bool try_pop(K& key) {
    sl_node_t *first, *next;
    bool result;

    first = head;

    while(1) {
      do {
        first = unset_mark(first->next[0]);
        next = first->next[0];
      } while(next && is_marked(next));

      if (next && !ATOMIC_CAS_MB(&first->next[0], next, set_mark(next))) {
      } else {
        break;
      }
    }

    result = (first->next[0] != NULL);
    if (!result) {
      return 0;
    }

    key = (first->key);
    mark_node_ptrs(first);

    fraser_search(key, NULL, NULL, first);
    sl_delete_node(first);

    return result;
  }

  K try_pop(void) {
    sl_node_t *first, *next;
    bool result;

    first = head;

    while(1) {
      do {
        first = unset_mark(first->next[0]);
        next = first->next[0];
      } while(next && is_marked(next));

      if (next && !ATOMIC_CAS_MB(&first->next[0], next, set_mark(next))) {
      } else {
        break;
      }
    }

    result = (first->next[0] != NULL);
    if (!result) {
      return first->key;
    }

    K key = (first->key);
    mark_node_ptrs(first);

    fraser_search(key, NULL, NULL, first);
    sl_delete_node(first);

    return key;
  }

  // SCANHEIGHT is what height to start spray at; must be >= 0
  #define SCANHEIGHT floor_log_2(n)+1
  // SCANMAX is scanlength at the top level; must be > 0
  #define SCANMAX floor_log_2(n)+1
  // SCANINC is the amount to increase scan length at each step; can be any integer
  #define SCANINC 0
  //SCANSKIP is # of levels to go down at each step; must be > 0
  #define SCANSKIP 1

  bool try_pop_spray(K& key, unsigned int n, sl_node_t **removed) {
    sl_node_t *cur;

retry:

    while (1) {
      sl_node_t *next;
      int scanlen;
      int height = SCANHEIGHT;
      int scanmax = SCANMAX;
      int scan_inc = SCANINC;
      int i = height;
      int dummy = 0;

      cur = head;

      while(1) {
        int r = _MarsagliaXOR();
        scanlen = r % (scanmax+1);

        while (dummy < n*floor_log_2(n)/2 && scanlen > 0) {
          dummy += (1 << i);
          scanlen--;
        }

        while (scanlen > 0 && cur->next[i]) { // Step right //here: cur->next[0], or cur->next[i]??
          sl_node_t *left = cur, *left_next = cur->next[i];
          if (is_marked(left_next)) goto retry;

          sl_node_t *right = left_next;
          while (1) {
              sl_node_t *right_next = right->next[i];
              if (!is_marked(right_next))
                break;
              right = unset_mark(right_next);
          }
          if (left_next != right) {
            if (!ATOMIC_CAS_MB(&left->next[i], left_next, right)) goto retry;
            for (sl_node_t *t = left_next; t != right; t = unset_mark(t->next[i]))
              t->next[i] = set_dead(t->next[i]);
          }
          cur = right;
          scanlen--;
        }

        // Got to end of list, maybe it's empty and maybe not.  Try a
        // normal pop() to be sure.
        if (!cur->next[0]) {
          return try_pop(key);
        }

        scanmax += scan_inc;

        if (i == 0) break;
        if (i <= SCANSKIP) { i = 0; continue; } // need to guarantee bottom level gets scanned
        i -= SCANSKIP;
      }

      if (cur == head) // still in dummy range
        return false; // TODO: clean instead? something else?

      for (next = cur->next[0]; is_marked(next) && next; ) {
        cur = unset_mark(next); // Find first non-deleted node
        next = cur->next[0];
      }

      if (!next) return try_pop(key);

      if (ATOMIC_CAS_MB(&cur->next[0], next, set_mark(next)))
        break;
    }

    key = (cur->key);
    mark_node_ptrs(cur);

    // Store nodes in local list for later reclamation
    cur->dummy = *removed;
    *removed = cur;

    return true;
  }
};

template<class Comparer, typename K>
Runtime::MM::ListNodeHeap LockFreeSkipList<Comparer,K>::heap[3];

template<class Comparer, typename K>
__thread unsigned long LockFreeSkipList<Comparer,K>::seeds[3];

template<class Comparer, typename K>
__thread bool LockFreeSkipList<Comparer,K>::seeds_init;

template<class Comparer, typename K>
__thread int LockFreeSkipList<Comparer,K>::spray_seed;

template<class Comparer, typename K>
__thread bool LockFreeSkipList<Comparer,K>::spray_seed_init;

template<class Comparer, typename K>
class SprayList : public LockFreeSkipList<Comparer, K> {

  typedef SkipListNode<K> sl_node_t;

  Runtime::PerThreadStorage<sl_node_t*> removedNodes;

  static bool node_linked(sl_node_t *n) {
    for (int i = n->toplevel - 1; i >= 0; i--) {
      if (!LockFreeSkipList<Comparer,K>::is_dead(n->next[i]))
        return true;
    }
    return false;
  }

  // check if nodes removed by spray have become unlinked
  // in the mean time, and reclaim them if so
  void cleanup(sl_node_t **head) {
    sl_node_t **prev = head;
    sl_node_t *n = *prev;
    sl_node_t *s = NULL;
    int maxlev = 0;

    while (n) {
      sl_node_t *next = n->dummy;

      if (!node_linked(n)) {
        *prev = n->dummy;
        this->sl_delete_node(n);
      } else {
        prev = &n->dummy;
        if (n->toplevel > maxlev) {
          maxlev = n->toplevel;
          s = n;
        }
      }
      n = next;
    }
    if (s)
      this->fraser_search(s->key, NULL, NULL, s);
  }

public:

  bool try_pop(K& key) {
    unsigned int n = Galois::getActiveThreads();
    sl_node_t **removed = removedNodes.getLocal();

    int r = LockFreeSkipList<Comparer,K>::_MarsagliaXOR();
    if (n == 1 || (r % n) == 0) { // n == 1 is equivalent to Lotan-Shavit delete_min
      cleanup(removed);
      return LockFreeSkipList<Comparer,K>::try_pop(key);
    }

    return this->try_pop_spray(key, n, removed);
  }
};


// MultiQueue, by Hamza Rihani, Peter Sanders, Roman Dementiev
// http://arxiv.org/abs/1411.1209
template<class Comparer, typename K, int c>
class MultiQueue {

private:
  LockFreeSkipList<Comparer, K> *Q;
  Comparer compare;
  int nQ;

public:
  MultiQueue() : nQ(Galois::getActiveThreads() * c) {
    Q = new LockFreeSkipList<Comparer, K>[nQ];
  }

  ~MultiQueue() {
    delete[] Q;
  }

  bool push(const K& key) {
    int q = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
    return Q[q].push(key);
  }

  bool try_pop(K& key) {

    while (true) {
      int q0 = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
      int q1 = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;

      if (q0 == q1) continue;

      SkipListNode<K> *first0 = Q[q0].peek_pop();
      SkipListNode<K> *first1 = Q[q1].peek_pop();
      bool gotit;

      if (!first0 && !first1) {
        for (int i = 0; i < nQ; i++) {
          if (Q[i].try_pop(key)) return true;
        }
        return false;
      } else if (!first0) {
        gotit = Q[q1].complete_pop(first1, key);
      } else if (!first1) {
        gotit = Q[q0].complete_pop(first0, key);
      } else if (compare(first0->key, first1->key)) {
        gotit = Q[q1].complete_pop(first1, key);
      } else {
        gotit = Q[q0].complete_pop(first0, key);
      }

      if (gotit) return true;
    }
  }
};

// MultiQueue, by Hamza Rihani, Peter Sanders, Roman Dementiev
// http://arxiv.org/abs/1411.1209
template<class Comparer, typename K, int c>
class HeapMultiQueue {

private:
  typedef boost::heap::d_ary_heap<K, boost::heap::arity<8>, boost::heap::compare<Comparer>> DAryHeap;
  struct Heap {
     Runtime::LL::SimpleLock<true> lock;
     K min;
     DAryHeap heap;
  };
  Runtime::LL::CacheLineStorage<Heap> *Q;
  Comparer compare;
  int nQ;
  K emptyK;

public:
  HeapMultiQueue() : nQ(Galois::getActiveThreads() * c) {
    Q = new Runtime::LL::CacheLineStorage<Heap>[nQ];
    memset(reinterpret_cast<void*>(&emptyK), 0xff, sizeof(emptyK));
    for (int i = 0; i < nQ; i++) {
      Q[i].data.min = emptyK;
      Q[i].data.heap.emplace(Q[i].data.min);
    }
  }

  bool push(const K& key) {
    Heap* h;
    int i;

    do {
      i = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
      h = &Q[i].data;
    } while (!h->lock.try_lock());

    h->heap.emplace(key);
    h->min = h->heap.top();
    h->lock.unlock();
    return true;
  }

  bool try_pop(K& key) {
    Heap *hi, *hj;
    int i, j;

    do {
      i = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
      hi = &Q[i].data;

      j = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
      hj = &Q[j].data;

      if (i == j) continue;

      if (compare(hi->min, hj->min))
        hi = hj;
    } while (!hi->lock.try_lock());

    if (hi->heap.size() == 1) {
      hi->lock.unlock();
      for (j = 1; j < nQ; j++) {
        hi = &Q[(i + j) % nQ].data;
        if (hi->min == emptyK) continue;
        hi->lock.lock();
        if (hi->heap.size() > 1)
          goto deq;
        hi->lock.unlock();
      }
      // empty
      return false;
    }

deq:
    key = hi->heap.top();
    hi->heap.pop();
    hi->min = hi->heap.top();
    hi->lock.unlock();
    return true;
  }
};

template<class Comparer, typename K, bool perPackage>
class DistQueue {

private:
  LockFreeSkipList<Comparer, K> *Q;
  Comparer compare;
  int nQ;

public:
  DistQueue() : nQ(Galois::getActiveThreads()) {
    unsigned nThr = Galois::getActiveThreads();

    nQ = perPackage ? Galois::Runtime::LL::getMaxPackageForThread(nThr - 1) + 1 : nThr;
    Q = new LockFreeSkipList<Comparer, K>[nQ];
  }

  ~DistQueue() {
    delete[] Q;
  }

  bool push(const K& key) {
    unsigned tid = Galois::Runtime::LL::getTID();
    unsigned qid = perPackage ? Galois::Runtime::LL::getPackageForThread(tid) : tid;

    return Q[qid].push(key);
  }

  bool try_pop(K& key) {
    unsigned tid = Galois::Runtime::LL::getTID();
    unsigned qid = perPackage ? Galois::Runtime::LL::getPackageForThread(tid) : tid;

    while (true) {
      SkipListNode<K> *min_node = 0;
      int q;

      for (int i = 0; i < nQ; i++) {
        int curq = (qid + i) % nQ;
        SkipListNode<K> *n = Q[curq].peek_pop();

        if (!n) continue;
        if (!min_node || compare(min_node->key, n->key)) {
          min_node = n;
          q = curq;
        }
      }

      if (!min_node)
        return false;
      if (Q[q].complete_pop(min_node, key))
        return true;
    }
  }
};

template<class Comparer, typename K, bool prioSteal>
class LocalPQ {

private:
  LockFreeSkipList<Comparer, K> *Q;
  Runtime::PerThreadStorage<unsigned int> current;
  Comparer compare;
  int nQ;

public:
  LocalPQ() : nQ(Galois::getActiveThreads()) {
    Q = new LockFreeSkipList<Comparer, K>[nQ];
  }

  ~LocalPQ() {
    delete[] Q;
  }

  bool push(const K& key) {
    unsigned tid = Galois::Runtime::LL::getTID();
    return Q[tid].push(key);
  }

  bool try_pop(K& key) {
    unsigned tid = Galois::Runtime::LL::getTID();
    int n = (*current.getLocal() += 1);

    if ((n & (1<<10)) != 0) {
        if (Q[tid].try_pop(key))
            return true;
    }

    while (prioSteal) {
        int q0 = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
        int q1 = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;

        if (q0 == q1) continue;

        SkipListNode<K> *first0 = Q[q0].peek_pop();
        SkipListNode<K> *first1 = Q[q1].peek_pop();
        bool gotit;

        if (!first0 && !first1) {
            break;
        } else if (!first0) {
            gotit = Q[q1].complete_pop(first1, key);
        } else if (!first1) {
            gotit = Q[q0].complete_pop(first0, key);
        } else if (compare(first0->key, first1->key)) {
            gotit = Q[q1].complete_pop(first1, key);
        } else {
            gotit = Q[q0].complete_pop(first0, key);
        }
        if (gotit)
            return true;
    }

    tid = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
    for (int i = 0; i < nQ; i++) {
      if (Q[(tid + i) % nQ].try_pop(key))
        return true;
    }

    return false;
  }
};

template<class Comparer, typename K>
class SwarmPQ {

private:
  LockFreeSkipList<Comparer, K> *Q;
  int nQ;

public:
  SwarmPQ() : nQ(Galois::getActiveThreads()) {
    Q = new LockFreeSkipList<Comparer, K>[nQ];
  }

  ~SwarmPQ() {
    delete[] Q;
  }

  bool push(const K& key) {
    unsigned tid = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
    return Q[tid].push(key);
  }

  bool try_pop(K& key) {
    unsigned tid = Galois::Runtime::LL::getTID();

    return Q[tid].try_pop(key);
  }
};

template<class Comparer, typename K>
class HeapSwarmPQ {

private:
  typedef boost::heap::d_ary_heap<K, boost::heap::arity<8>, boost::heap::compare<Comparer>> DAryHeap;
  struct Heap {
     Runtime::LL::SimpleLock<true> lock;
     DAryHeap heap;
  };
  Runtime::LL::CacheLineStorage<Heap> *Q;
  int nQ;

public:
  HeapSwarmPQ() : nQ(Galois::getActiveThreads()) {
    Q = new Runtime::LL::CacheLineStorage<Heap>[nQ];
  }

  bool push(const K& key) {
    unsigned tid;

    while (1) {
      tid = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
      if (Q[tid].data.lock.try_lock()) break;
    }
    Q[tid].data.heap.emplace(key);
    Q[tid].data.lock.unlock();
    return true;
  }

  bool try_pop(K& key) {
    unsigned tid = Galois::Runtime::LL::getTID();

    Q[tid].data.lock.lock();
    if (!Q[tid].data.heap.empty()) {
        key = Q[tid].data.heap.top();
        Q[tid].data.heap.pop();
        Q[tid].data.lock.unlock();
        return true;
    }
    Q[tid].data.lock.unlock();

    return false;
  }
};


template<class Comparer, class Hasher, typename K>
class PartitionPQ {

private:
  LockFreeSkipList<Comparer, K> *Q;
  Hasher hash;
  int nQ;

public:
  PartitionPQ() : nQ(Galois::getActiveThreads()) {
    Q = new LockFreeSkipList<Comparer, K>[nQ];
  }

  ~PartitionPQ() {
    delete[] Q;
  }

  bool push(const K& key) {
    static const unsigned long s = 2654435769ull ;
    const unsigned long h = hash(key);
    const unsigned long q = (h * s) & 0xffffffff ;
    return Q[q % nQ].push(key);
  }

  bool try_pop(K& key) {
    unsigned tid = Galois::Runtime::LL::getTID();

    if (Q[tid].try_pop(key))
        return true;

    tid = LockFreeSkipList<Comparer, K>::rand_range(nQ) - 1;
    for (int i = 0; i < nQ; i++) {
      if (Q[(tid + i) % nQ].try_pop(key))
        return true;
    }

    return false;
  }
};

// check cache alignment
template<typename K,typename V>
struct SkipListSetNode {
public:
  // dummy field which is used by the heap when the node is freed.
  // (without it, freeing a node would corrupt a field, possibly affecting
  // a concurrent traversal.)
  void* dummy;

  K key;
  V val;

  int toplevel;
  SkipListSetNode* volatile next[0];

  void init(int _t, SkipListSetNode *_next) {
    toplevel = _t;
    for (int i = 0; i < _t; i++)
      next[i] = _next;
  }

};

template<class Comparer, typename K, typename V>
class LockFreeSkipListSet {

private:
  typedef SkipListSetNode<K,V> sl_node_t;

  static Runtime::MM::ListNodeHeap heap[3];
  Runtime::TerminationDetection& term;
  Comparer compare;

public:
  sl_node_t* head;
  uint8_t levelmax;

  static inline bool is_marked(sl_node_t* i)
  {
    return ((uintptr_t)i & (uintptr_t)0x01) != 0;
  }

  static inline sl_node_t* unset_mark(sl_node_t* i)
  {
    return (sl_node_t *)((uintptr_t)i & ~(uintptr_t)0x03);
  }

  static inline sl_node_t* set_mark(sl_node_t* i)
  {
    return (sl_node_t *)((uintptr_t)i | (uintptr_t)0x01);
  }

  //Marsaglia's xorshf generator
  static inline unsigned long xorshf96(unsigned long* x, unsigned long* y, unsigned long* z)  //period 2^96-1
  {
    unsigned long t;
    (*x) ^= (*x) << 16;
    (*x) ^= (*x) >> 5;
    (*x) ^= (*x) << 1;

    t = *x;
    (*x) = *y;
    (*y) = *z;
    (*z) = t ^ (*x) ^ (*y);

    return *z;
  }

  static __thread unsigned long seeds[3];
  static __thread bool seeds_init;

public:
  static inline long rand_range(long r)
  {
    if (!seeds_init) {
      int fd = open("/dev/urandom", O_RDONLY);
      if (read(fd, seeds, 3 * sizeof(unsigned long)) < 0) {
        perror("read");
        exit(1);
      }
      close(fd);
      seeds_init = true;
    }
    long v = xorshf96(seeds, seeds + 1, seeds + 2) % r;
    v++;
    return v;
  }

private:
  void mark_node_ptrs(sl_node_t *n)
  {
    sl_node_t *n_next;
    int i;

    for (i=n->toplevel-1; i>=0; i--)
    {
      do
      {
        n_next = n->next[i];
        if (is_marked(n_next))
        {
          break;
        }
      } while (!ATOMIC_CAS_MB(&n->next[i], n_next, set_mark(n_next)));
    }
  }

  static int floor_log_2(unsigned int n) {
    int pos = 0;
    if (n >= 1<<16) { n >>= 16; pos += 16; }
    if (n >= 1<< 8) { n >>=  8; pos +=  8; }
    if (n >= 1<< 4) { n >>=  4; pos +=  4; }
    if (n >= 1<< 2) { n >>=  2; pos +=  2; }
    if (n >= 1<< 1) {           pos +=  1; }
    return ((n == 0) ? (-1) : pos);
  }

  inline sl_node_t *sl_new_node(int levelmax, sl_node_t *next) {
    int e = term.getEpoch() % 3;
    sl_node_t *node = reinterpret_cast<sl_node_t *>(heap[e].allocate(sizeof(sl_node_t) + levelmax*sizeof(sl_node_t*), levelmax-1));
    node->init(levelmax, next);
    return node;
  }

  inline sl_node_t *sl_new_node_key(K key, V val, int levelmax) {
    sl_node_t *node = sl_new_node(levelmax, 0);
    node->key = key;
    node->val = val;
    return node;
  }

  inline void sl_delete_node(sl_node_t *n) {
    int e = (term.getEpoch() + 2) % 3;
    heap[e].deallocate(n, n->toplevel-1);
  }

public:

  LockFreeSkipListSet() : term(Runtime::getSystemTermination()), levelmax(23) {
    sl_node_t *min, *max;

    max = sl_new_node(levelmax, NULL);
    min = sl_new_node(levelmax, max);

    head = min;
  }

  bool empty() const
  {
    return head->next[0]->next[0] == 0;
  }

  void fraser_search(K key, sl_node_t **left_list, sl_node_t **right_list, sl_node_t *dead)
  {
    sl_node_t *left, *left_next, *right, *right_next;
    int i;

  retry:
    left = head;
    for (i = (dead ? dead->toplevel : levelmax) - 1; i >= 0; i--)
    {
      left_next = left->next[i];
      if (is_marked(left_next))
        goto retry;
      /* Find unmarked node pair at this level */
      for (right = left_next; ; right = right_next)
      {
        /* Skip a sequence of marked nodes */
        while(1)
        {
          right_next = right->next[i];
          if (!is_marked(right_next))
            break;
          right = unset_mark(right_next);
        }
        /* Ensure left and right nodes are adjacent */
        if ((left_next != right) &&
            !ATOMIC_CAS_MB(&left->next[i], left_next, right))
          goto retry;
        /* When deleting, we have to keep going until we find our target node, or until
           we observe it has been deleted (right->key > key).  Once this happens, however,
           we need to descend to a node whose key is smaller than our target's, otherwise
           we might miss our target in the level below.  (Consider N1 and N2 with the same
           key, where on level 1, N1 -> N2 but on level 0, N2 -> N1; if when looking for N2
           we descend at N1, we'll miss N2 at level 0.) */
        if (!right_next || !compare(key, right->key))
          break;
        left = right;
        left_next = right_next;
      }
      if (left_list != NULL)
        left_list[i] = left;
      if (right_list != NULL)
        right_list[i] = right;
    }
  }

  int get_rand_level()
  {
    int i, level = 1;
    for (i = 0; i < levelmax - 1; i++)
      {
        if ((rand_range(100)-1) < 50)
          level++;
        else
          break;
      }
    /* 1 <= level <= *levelmax */
    return level;
  }

  V get(const K& key) {
    sl_node_t *succs[levelmax], *preds[levelmax];

    fraser_search(key, preds, succs, NULL);
    if (succs[0]->next[0] && succs[0]->key == key)
      return succs[0]->val;
    return static_cast<V>(0);
  }

  sl_node_t* lower_bound(const K& key) {
    sl_node_t *succs[levelmax], *preds[levelmax];

    fraser_search(key, preds, succs, NULL);
    return succs[0];
  }

  bool pop(sl_node_t* node) {
    sl_node_t *first, *next;
    bool result;

    next = node->next[0];
    if (!next || !ATOMIC_CAS_MB(&node->next[0], next, set_mark(next)))
      return false;

    mark_node_ptrs(node);

    fraser_search(node->key, NULL, NULL, node);
    sl_delete_node(node);

    return true;
  }

  bool push(const K& key, const V& val)
  {
    sl_node_t *newn, *new_next, *pred, *succ, *succs[levelmax], *preds[levelmax];
    int i, result = 0;

    newn = sl_new_node_key(key, val, get_rand_level());

  retry:
    fraser_search(key, preds, succs, NULL);
    if (succs[0]->next[0] && succs[0]->key == key)
    {                             /* Value already in list */
      result = 0;
      sl_delete_node(newn);
      goto end;
    }

    for (i = 0; i < newn->toplevel; i++)
    {
      newn->next[i] = succs[i];
    }

    /* Node is visible once inserted at lowest level */
    if (!ATOMIC_CAS_MB(&preds[0]->next[0], succs[0], newn))
    {
      goto retry;
    }

    for (i = 1; i < newn->toplevel; i++)
    {
      while (1)
      {
        pred = preds[i];
        succ = succs[i];
        new_next = newn->next[i];
        /* Give up if pointer is marked */
        if (is_marked(new_next))
          goto success;
        /* Update the forward pointer if it is stale, which can happen
           if we called search again to update preds and succs. */
        if (new_next != succ && !ATOMIC_CAS_MB(&newn->next[i], new_next, succ))
          goto success;
        /* We retry the search if the CAS fails */
        if (ATOMIC_CAS_MB(&pred->next[i], succ, newn)) {
          if (is_marked(newn->next[i])) {
            fraser_search(key, NULL, NULL, newn);
            goto success;
          }
          break;
        }

        fraser_search(key, preds, succs, NULL);
      }
    }

  success:
    result = 1;

  end:
    return result;
  }

};

template<class Comparer, typename K, typename V>
Runtime::MM::ListNodeHeap LockFreeSkipListSet<Comparer,K,V>::heap[3];

template<class Comparer, typename K, typename V>
__thread unsigned long LockFreeSkipListSet<Comparer,K,V>::seeds[3];

template<class Comparer, typename K, typename V>
__thread bool LockFreeSkipListSet<Comparer,K,V>::seeds_init;


template<typename K, class Indexer, int Rlx>
class kLSMQ {
    // kpq::k_lsm key MUST be unsigned
    kpq::k_lsm<unsigned long, K, Rlx> pq;
    Indexer indexer;

public:
    bool push(const K& key) {
      pq.insert((unsigned long)indexer(key), key);
      return true;
    }

    bool try_pop(K& key) {
      return pq.delete_min(key);
    }
};

}
} // end namespace Galois

#endif
