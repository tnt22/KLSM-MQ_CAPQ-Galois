/*
 *  Copyright 2014 Jakob Gruber
 *
 *  This file is part of kpqueue.
 *
 *  kpqueue is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  kpqueue is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with kpqueue.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __GENERIC_DIST_LSM_LOCAL_H
#define __GENERIC_DIST_LSM_LOCAL_H

#include <atomic>

#include "block_storage.h"
#include "item.h"
#include "counters.h"
#include "mm.h"
#include "xorshf96.h"

namespace kpq
{

template <class K, class V, int Rlx, class PQ>
class generic_dist_lsm;

template <class K, class V, int Rlx, class PQ>
class generic_dist_lsm_local
{
public:
    generic_dist_lsm_local();
    virtual ~generic_dist_lsm_local();

    void insert(const K &key,
                const V &val,
                PQ *slsm);
    bool delete_min(generic_dist_lsm<K, V, Rlx, PQ> *parent,
                    V &val);
    bool delete_min(generic_dist_lsm<K, V, Rlx, PQ> *parent,
                    K &key, V &val);
    /** Iterates through local items and returns the best one found.
     *  In the process of finding the minimal item, unowned items
     *  in each block are removed and block merges are performed if possible.
     *  Used internally by delete_min() and by the k-lsm's delete_min()
     *  operation. */
    void peek(typename block<K, V>::peek_t &best);

    /** Performs a peek without mutating blocks May be called from other threads. */
    void safe_peek(typename block<K, V>::peek_t &best);

    /** Attempts to copy items from a random other thread's local clsm,
     *  and returns the number of items copied. */
    int spy(class generic_dist_lsm<K, V, Rlx, PQ> *parent);
    int spy(generic_dist_lsm_local<K, V, Rlx, PQ> *victim);

    bool empty() const { return m_head.load(std::memory_order_relaxed) == nullptr; }

    void print() const;

private:
    /** The internal insertion, used both in the public insert() and in spy(). */
    void insert(item<K, V> *it,
                const version_t version,
                PQ *slsm);

    /**
     * Inserts new_block into the linked list of blocks, merging with
     * same size blocks until no two blocks in the list have the same size.
     */
    void merge_insert(block<K, V> *const new_block,
                      PQ *slsm);

private:
    std::atomic<block<K, V> *> m_head; /**< The largest  block. */
    block<K, V>               *m_tail; /**< The smallest block. */
    block<K, V>               *m_spied;

    block_storage<K, V, 4> m_block_storage;
    item_allocator<item<K, V>, typename item<K, V>::reuse> m_item_allocator;

    /** Caches the previously peeked item in case we can short-circuit and simply
     *  return it. */
    typename block<K, V>::peek_t m_cached_best;

    xorshf96 m_gen;
};

#include "generic_dist_lsm_local_inl.h"

}

#endif /* __GENERIC_DIST_LSM_LOCAL_H */
