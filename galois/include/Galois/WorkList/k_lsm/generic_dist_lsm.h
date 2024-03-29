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

#ifndef __GENERIC_DIST_LSM_H
#define __GENERIC_DIST_LSM_H

#include "generic_dist_lsm_local.h"

namespace kpq
{

template <class K, class V, int Rlx, class PQ>
class generic_dist_lsm
{
    friend int generic_dist_lsm_local<K, V, Rlx, PQ>::spy(generic_dist_lsm<K, V, Rlx, PQ> *parent);

public:

    /**
     * Inserts a new item into the local LSM.
     */
    void insert(const K &key);
    void insert(const K &key,
                const V &val);

    /**
     * A special version of insert for use by the k-lsm. Acts like a standard
     * insert until the largest block exceeds the relaxation size limit, at which
     * point the block is inserted into the shared lsm instead.
     */
    void insert(const K &key,
                const V &val,
                PQ *slsm);

    /**
     * Attempts to remove the locally (i.e. on the current thread) minimal item.
     * If the local LSM is empty, we try to copy items from another active thread.
     * In case the local LSM is still empty, false is returned.
     * If a locally minimal element is successfully found and removed, true is returned.
     */
    bool delete_min(V &val);
    bool delete_min(K &key, V &val);
    void find_min(typename block<K, V>::peek_t &best);

    int spy();

    void print();

    void init_thread(const size_t) const { }
    constexpr static bool supports_concurrency() { return true; }

private:
    thread_local_ptr<generic_dist_lsm_local<K, V, Rlx, PQ>> m_local;
};

#include "generic_dist_lsm_inl.h"

}

#endif /* __GENERIC_DIST_LSM_H */
