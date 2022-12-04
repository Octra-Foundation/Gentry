/*





MIT License

Copyright (c) 2022 Octra Foundation

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.




*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <math.h>
#include <stdexcept>
#include <algorithm>
#include <emmintrin.h>

#include "bit_hash_map.h"

typedef unsigned char  UC;
typedef unsigned short US;
typedef unsigned int   UI;

using namespace std;

template <class Template> void alloc(Template*& next_block, int block) {
    next_block = (Template*)calloc(block, sizeof(Template));
    if (!next_block) {
        static_cast<void>(fprintf(stderr, "out of memory\n")), exit(1);
    }
}

std::ifstream::pos_type filesize(const char* f_input) {
    std::ifstream in(f_input, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

void quit(const char* message = 0) {
  throw message;
}

int equals(const char* string_a, const char* string_b) {
  assert(string_a && string_b);
  while (*string_a && *string_b) {
    if (std::tolower(*string_a) != std::tolower(*string_b)) return 0;
    ++string_a;
    ++string_b;
  }
  return *string_a ==* string_b;
}



#define GENERAL_TABLE_VALUE(x) (32768 / (x + x + 3))
#define BLOCK_VALUE(x) ((x & 1) * 2 + (x & 2) + (x >> 2 & 1) + (x >> 3 & 1) + (x >> 4 & 1) + (x >> 5 & 1) + (x >> 6 & 1) + (x >> 7 & 1) + 3)



class State_Machine {
public:
    State_Machine(int block);

    int next_block(int state) {
        assert(state >= 0 && state < state_number);
        return state_count[previous_state = state] >> 16;
    }

    void update(int rm_vector, int limit = 255) {
        int block = state_count[previous_state] & 255, next_block = state_count[previous_state] >> 14;
        if (block < limit) {
            ++state_count[previous_state];
            state_count[previous_state] += ((rm_vector << 18) - next_block) * general_table[block] & 0xffffff00;
        }
    }

private:
    int state_number;
    int previous_state;
    int state_count[];
    static const int general_table[256];
};

const int State_Machine::general_table[256] = {
    GENERAL_TABLE_VALUE(0),
    GENERAL_TABLE_VALUE(1),
    GENERAL_TABLE_VALUE(2),
    GENERAL_TABLE_VALUE(253),
    GENERAL_TABLE_VALUE(254),
    GENERAL_TABLE_VALUE(255)
};

State_Machine::State_Machine(int block) : state_number(block), previous_state(0) {
    alloc(state_count, state_number);
    for (int i = 0; i < state_number; ++i) {
        int block = BLOCK_VALUE(i);
        state_count[i] = block << 28 | 6;
    }
}

State_Machine::State_Machine(int block) {
  state_number = block;
  alloc(state_count, block);
  memset(state_count, 0, block * sizeof(int));
}

#undef GENERAL_TABLE_VALUE
#undef BLOCK_VALUE



class Determinant {
public:
    Determinant();
    std::uint16_t next_block() {
        return directed_graph.next_block(previous_state << 8 | state[previous_state]);
    }

    void update(std::uint16_t rm_vector) {
        directed_graph.update(rm_vector, 90);
        previous_state = (previous_state << 8 | state[previous_state]) + rm_vector;
        state[previous_state >> 8] = previous_state & 255;
        previous_state >>= 8;
    }

private:
    std::uint16_t previous_state;
    State_Machine directed_graph;
    std::array<int, 256> state;
};

Determinant::Determinant() : previous_state(0), directed_graph(0x10000), state({}) {
    std::fill(state.begin(), state.end(), 0x66);
}


typedef enum { 
    ENCODE, 
    DECODE 
    } 
Mode;

class Encoder {
private:
    Determinant determinant;
    const Mode mod;
    FILE* __data;
    UI vector_x, vector_y;
    UI lambda;
public:
    Encoder(Mode io_mode, FILE* __str); // for ease of testing will be used with a data file
    void encode(int rm_vector);
    int decode();
    void alignment();

};

Encoder::Encoder(Mode io_mode, FILE* __str) : determinant(), mod(io_mode), __data(__str), vector_x(0),
vector_y(0xffffffff), lambda(0) {
    if (mod == DECODE) {
        for (int i = 0; i < 4; ++i) {
            int byte = getc(__data);
            if (byte == EOF) byte = 0;
            lambda = (lambda << 8) + (byte & 0xff);
        }
    }
}

inline void Encoder::encode(int rm_vector) {
    const UI next_block = determinant.next_block();
    assert(next_block <= 0xffff);
    assert(rm_vector == 0 || rm_vector == 1);
    const UI y_minus_x = vector_y - vector_x;
    const UI symmetric = vector_x + (y_minus_x >> 16) * next_block + ((y_minus_x & 0xffff) * next_block >> 16);
    assert(symmetric >= vector_x && symmetric < vector_y);
    if (rm_vector) {
        vector_y = symmetric;
    }
    else {
        vector_x = symmetric + 1;
        determinant.update(rm_vector);
        while (((vector_x ^ vector_y) & 0xff000000) == 0) {
            putc(vector_y >> 24, __data);
            vector_x <<= 8;
            vector_y = (vector_y << 8) + 255;
        }
    }
}


inline int Encoder::decode() {
    const UI next_block = determinant.next_block();
    assert(next_block <= 0xffff);
    const UI y_minus_x = vector_y - vector_x;
    const UI symmetric = vector_x + (y_minus_x >> 16) * next_block + ((y_minus_x & 0xffff) * next_block >> 16);
    assert(symmetric >= vector_x && symmetric < vector_y);
    int rm_vector = 0;
    if (lambda <= symmetric) {
        rm_vector = 1;
        vector_y = symmetric;
    }
    else {
        vector_x = symmetric + 1;
    }
    determinant.update(rm_vector);

    while (((vector_x ^ vector_y) & 0xff000000) == 0) {
        vector_x <<= 8;
        vector_y = (vector_y << 8) + 255;
        int byte = getc(__data);
        if (byte == EOF) byte = 0;
        lambda = (lambda << 8) + byte;
    }
    return rm_vector;
}


void Encoder::alignment() {
    if (mod == ENCODE) {
        const UI x = vector_x;
        const UI y = vector_y;
        while (((x ^ y) & 0xff000000) == 0) {
            putc(y >> 24, __data);
            vector_x <<= 8;
            vector_y = (y << 8) + 255;
        }
        putc(y >> 24, __data);
    }
}

class String: public Array<char> {
public:
  const char* c_str() const {
    return &(*this)[0];
  }
  
  void operator = (const char* s) {
    resize(strlen(s) + 1);
    strcpy(&(*this)[0], s);
  }
  
  void operator += (const char* s) {
    pop_back();
    while (*s) push_back(*s++);
    push_back(0);
  }
  String(const char* s = ""): Array<char>(1) {
    (*this) += s;
  }
};

template <class T, int MOVE = 0> class Array {
private:
  int   local_size;
  int   __res;
  char  *memory;
  T*    f_data; 
  void create(int i);

public:
  explicit Array(int i = 0) {
    create(i);
  }

  ~Array();

  T& operator[](int i) {
#ifndef NDEBUG
    if (i < 0 || i >= local_size) quit();
#endif
    return f_data[i];
  }

  const T& operator[](int i) const {
#ifndef NDEBUG
    if (i < 0 || i >= local_size) quit();
#endif
    return f_data[i];
  }

  int size() const {
    return local_size;
  }

  void resize(int i);

  void pop_back() {
    if (local_size > 0) 
      --local_size;
  }

  void push_back(const T& x);
private:
  Array(const Array&);
  Array& operator=(const Array&);
};

template<class T, int MOVE> void Array<T, MOVE>::resize(int i) {
  if (i <= __res) {
    local_size = i;
    return;
  }

  char *ptr_memory = memory;
  T *__st = f_data;
  int complex_size = local_size;
  create(i);
  if (ptr_memory) {
    if (__st) {
      memcpy(f_data, __st, sizeof(T) * min(i, complex_size));
    }
    free(ptr_memory);
  }
}

template<class T, int MOVE> void Array<T, MOVE>::create(int i) {
  local_size = __res = i;
  if (i <= 0) {
    f_data = 0;
    memory = 0;
    return;
  }

  const int __size = MOVE + local_size * sizeof(T);
  memory = (char*)calloc(__size, 1);
  if (!memory) quit();

  f_data = (MOVE ? (T*)(memory + MOVE - (((long)memory) & (MOVE - 1))) : (T*)memory);
  assert((char*)f_data >= memory && (char*)f_data <= memory + MOVE);
}

template<class T, int MOVE> Array<T, MOVE>::~Array() {
  free(memory);
}

template<class T, int MOVE> void Array<T, MOVE>::push_back(const T& x) {
  if (local_size == __res) {
    int complex_size = local_size;
    resize(max(1, local_size * 2));
    local_size = complex_size;
  }
  f_data[local_size++] = x;
}

class Random {
  Array<UI> table;
  int i;
    
public:
  Random(): table(64) {
    table[0] = 123456789;
    table[1] = 987654321;
      
    for(int j=  0; j < 62; j++) 
        table[j + 2] = table[j + 1] * 11 + table[j] * 23 / 16;
    i = 0;
  }
    
  UI operator()() {
    return ++i, table[i&63] = table[i - 24&63] ^ table[i - 55&63];
  }
    
} random;

class Ilog {
  Array<UC> t;
public:
  int operator()(U16 x) const {
      return t[x];
  }
  Ilog();
} ilog;

Ilog::Ilog() : t(65536) {
  uint32_t x = 14155776;
  for (int i = 2; i < 65536; ++i) {
    x += 774541002 / (i * 2 - 1);  // numerator is 2^29/ln 2
    t[i] = x >> 24;
  }
}

inline int llog(UI x) {
  if (x >= 0x1000000) return 256 + ilog(x >> 16);
  else if (x >= 0x10000) return 128 + ilog(x >> 8);
  else return ilog(x);
}

int scalar_product(short *R, short *Q, int n) {
    __m128i sum = _mm_set1_epi32(0);
    n = (n + 7) & -8;
    for (int i = 0; i < n; i += 8) {
        __m128i R_vec = _mm_loadu_si128((__m128i*)&R[i]);
        __m128i Q_vec = _mm_loadu_si128((__m128i*)&Q[i]);
        __m128i prod = _mm_mullo_epi16(R_vec, Q_vec);
        sum = _mm_add_epi32(sum, prod);
    }
    
    int result = _mm_extract_epi32(sum, 0) + _mm_extract_epi32(sum, 1) + _mm_extract_epi32(sum, 2) + _mm_extract_epi32(sum, 3);
    return result >> 8;
}


template <int Bit>
class BitHashMap {
public:
  static const int SearchLimit = 8;
  explicit BitHashMap(int size)
      : elements_(size * Bit), size_(size - 1) {
    assert(Bit >= 2 && size > 0 && (size & (size - 1)) == 0);
  }
  unsigned char* operator[](unsigned int index);

private:
  std::array<unsigned char, 64> elements_;
  unsigned int size_;
};

template <int Bit>
unsigned char* BitHashMap<Bit>::operator[](unsigned int index) {
  int check = (index >> 16 ^ index) & 0xffff;
  index = index * SearchLimit & size_;
  unsigned char* element;
  unsigned short* context;
  int j;
  for (j = 0; j < SearchLimit; ++j) {
    element = &elements_[(index + j) * Bit];
    context = reinterpret_cast<unsigned short*>(element);
    if (element[2] == 0) *context = check;
    if (*context == check) break;
  }
  if (j == 0) return element + 1;
  static unsigned char tmp[Bit];
  if (j == SearchLimit) {
    --j;
    memset(tmp, 0, Bit);
    *reinterpret_cast<unsigned short*>(tmp) = check;
    if (SearchLimit > 2 && elements_[(index + j) * Bit + 2] >
                               elements_[(index + j - 1) * Bit + 2])
      --j;
  } else
    memcpy(tmp, context, Bit);
  memmove(&elements_[(index + 1) * Bit], &elements_[index * Bit], j * Bit);
  memcpy(&elements_[index * Bit], tmp, Bit);
  return &elements_[index * Bit + 1];
}

int main() {
    State_Machine sm;
        for (int i = 0; i < sm.state_number; ++i) {
            sm.update(i, i % 256);
        }
    return 0;
}
