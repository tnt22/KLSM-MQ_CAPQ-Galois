#ifndef CAPQ_CLASS_H
#define CAPQ_CLASS_H

#include "Galois/WorkList/capq/capq_c_library/cppcapq.h"

namespace cpq
{

template <typename K, class Indexer>
class CA_PQ
{
private:
  CPPCAPQ<true,true,true> capq;
  Indexer indexer;

public:
  CA_PQ() {  }
  ~CA_PQ() {  }
  
  bool push(const K &key)
  {
    size_t ind = (size_t)indexer(key);
    K *keyCopy = new K();
    *keyCopy = key;
    capq.insert(ind, (size_t)keyCopy);
    return true;
  }

  bool try_pop(K &key)
  {
    size_t key_ret, val_ret;
    bool result = capq.delete_min(key_ret, val_ret);
    if(result) {
        key = *((K*)val_ret);
        delete ((K*)val_ret);
    }
    return result;
  }
};
} // namespace cpq
#endif
