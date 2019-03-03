/*
 *  Copyright 2015 Jakob Gruber & Ido Kessler
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

#ifndef __GENERIC_K_LSM_H
#define __GENERIC_K_LSM_H

#include "block.h"
#include "generic_dist_lsm.h"
#include "counters.h"

namespace kpq
{

/**
 * The k-lsm combines the distributed- and shared lsm data structures
 * in order to emphasize their respective strenghts. Items are initially
 * inserted into (thread-local) distributed lsm's until the relaxation
 * limit is reached, at which point the contained item's are inserted
 * into the shared lsm component.
 *
 * As always, K, V and Rlx denote, respectively, the key, value classes
 * and the relaxation parameter.
 * The PQ parameter denote the Priority Queue to use for the algorithm.
 * The PQ should implement void push(V& val) and bool try_pop(V& val);
 */

template <class K, class V, int Rlx, class Indexer, class PQ>
class generic_k_lsm
{
public:
  generic_k_lsm();
  virtual ~generic_k_lsm() {}

  void insert(const K &key);
  void insert(const K &key,
              const V &val);

  int delete_min(V &val);
  int delete_min(K &key, V &val);

  void init_thread(const size_t) const {}
  constexpr static bool supports_concurrency() { return true; }

private:
  Indexer indexer;
  class PQWraper
  {
  private:

  public:
    PQ pq;
  
    void insert(block<K, V> *block)
    {
      typename kpq::block<K, V>::peek_t item;
      V val;
      while ((item = block->peek()).m_item != nullptr)
        if (item.take(val))
          pq.push(val);
    }
  };
  generic_dist_lsm<K, V, Rlx, PQWraper> m_dist;
  PQWraper m_shared;
};

#include "generic_k_lsm_inl.h"

} // namespace kpq

#endif /* __GENERIC_K_LSM_H */
