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


#include "lambda.cpp" //pls use GMP (https://gmplib.org/)



// 8, 16, 32 bit unsigned types (adjust as appropriate)


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

class State_Machine {
protected:
    const int state_number;
    int previous_state;
    UI* state_count;
    static int general_table[256];
    
public:
    State_Machine(int block = 256);
    
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
};

int State_Machine::general_table[256] = { 0 };

State_Machine::State_Machine(int block) : state_number(block), previous_state(0) {
    alloc(state_count, state_number);
    
    for (int i = 0; i < state_number; ++i) {
        UI block = (i & 1) * 2 + (i & 2) + (i >> 2 & 1) + (i >> 3 & 1) + (i >> 4 & 1) + (i >> 5 & 1) + (i >> 6 & 1) + (i >> 7 & 1) + 3;
        state_count[i] = block << 28 | 6;
    }
    
    if (general_table[0] == 0) {
        for (int i = 0; i < 256; ++i) {
            general_table[i] = 32768 / (i + i + 3);
        }
    }
}

class Determinant {
    int previous_state;
    State_Machine directed_graph;
    int state[256];
public:
    Determinant();
    int next_block() {
        return directed_graph.next_block(previous_state << 8 | state[previous_state]);
    }

    void update(int rm_vector) {
        directed_graph.update(rm_vector, 90);
        int& state_numeration = state[previous_state];
        (state_numeration += state_numeration + rm_vector) &= 255;
        if ((previous_state += previous_state + rm_vector) >= 256)
            previous_state = 0;
    }
};

Determinant::Determinant() : previous_state(0), directed_graph(0x10000) {
    for (int i = 0; i < 0x100; ++i)
        state[i] = 0x66;
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
    const UI __symmetric = vector_x + ((vector_y - vector_x) >> 16) * next_block + ((vector_y - vector_x & 0xffff) * next_block >> 16);
    assert(__symmetric >= vector_x && __symmetric < vector_y);
    if (rm_vector) {
        vector_y = __symmetric;
    }
    else
        vector_x = __symmetric + 1;
        determinant.update(rm_vector);
        while (((vector_x ^ vector_y) & 0xff000000) == 0) {
            putc(vector_y >> 24, __data);
            vector_x <<= 8;
            vector_y = (vector_y << 8) + 255;
    }
}

inline int Encoder::decode() {
    const UI next_block = determinant.next_block();
    assert(next_block <= 0xffff);
    const UI __symmetric = vector_x + ((vector_y - vector_x) >> 16) * next_block + ((vector_y - vector_x & 0xffff) * next_block >> 16);
    assert(__symmetric >= vector_x && __symmetric < vector_y);
    int rm_vector = 0;
    if (lambda <= __symmetric) {
        rm_vector = 1;
        vector_y = __symmetric;
    }
    else
        vector_x = __symmetric + 1;
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
        while (((vector_x ^ vector_y) & 0xff000000) == 0) {
            putc(vector_y >> 24, __data);
            vector_x <<= 8;
            vector_y = (vector_y << 8) + 255;
        }
        putc(vector_y >> 24, __data);
    }
}



//////////////////////////// Random Generator ///////////////////////////////
// 32bit p-random number generator

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

///////////////////////////// ilog //////////////////////////////

// ilog(x) = round(log2(x) * 16), 0 <= x < 64K
class Ilog {
  Array<UC> t;
public:
  int operator()(U16 x) const {
      return t[x];
  }
  Ilog();
} ilog;

// Lookup-tab by num integration of 1/x
Ilog::Ilog(): t(65536) {
  UI x = 14155776;
  for (int i = 2; i < 65536; ++i) {
    x += 774541002 / (i * 2 - 1);  // numerator is 2^29/ln 2
    t[i] = x >> 24;
  }
}

// llog(x) accepts 32 bits
inline int llog(UI x) {
  if (x >= 0x1000000)
    return 256 + ilog(x >> 16);
  else if (x >= 0x10000)
    return 128 + ilog(x >> 8);
  else
    return ilog(x);
}



int main(int argc, char** argv) {
    /*
    EDT next
    */
    int byte;
    if (argv[1][0] == 'c') {
        Encoder coder(ENCODE, out);
        while ((byte = getc(in)) != EOF) {
           coder.encode(1);
            for (int i = 7; i >= 0; --i)
                coder.encode((byte >> i) & 1);
        }
        coder.encode(0);
        coder.alignment();
    }
    return 0;
}
