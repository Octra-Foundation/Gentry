#ifndef BIT_HASH_MAP_H
#define BIT_HASH_MAP_H

template <int _bit>
class bit_hash_map {
  enum {
  map = 8
  };

  Array<UC, 64> t;
  U32 size;
public:
  bit_hash_map(int i): t(i*_bit), size(i-1) {
    assert(_bit>=2 && i>0 && (i&(i-1))==0);
  }
  UC* operator[](U32 i);
};

#endif
