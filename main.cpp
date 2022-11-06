/*





MIT License

Copyright (c) [2022] [Octra Foundation]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files, to deal
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

typedef unsigned char  UC;
typedef unsigned short US;
typedef unsigned int   UI;
using namespace std;

template <class Template> void alloc(Template*& next_block, int block_count) {
    next_block = (Template*)calloc(block_count, sizeof(Template));
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
    State_Machine(int block_count = 256);
    
    int next_block(int state) {
        assert(state >= 0 && state < state_number);
        return state_count[previous_state = state] >> 16;
    }

    void update(int rm_vector, int limit = 255) {
        int block_count = state_count[previous_state] & 255, next_block = state_count[previous_state] >> 14;
        if (block_count < limit) {
            ++state_count[previous_state];
            state_count[previous_state] += ((rm_vector << 18) - next_block) * general_table[block_count] & 0xffffff00;
        }
    }
};

int State_Machine::general_table[256] = { 0 };

State_Machine::State_Machine(int block_count) : state_number(block_count), previous_state(0) {
    alloc(state_count, state_number);
    
    for (int i = 0; i < state_number; ++i) {
        UI block_count = (i & 1) * 2 + (i & 2) + (i >> 2 & 1) + (i >> 3 & 1) + (i >> 4 & 1) + (i >> 5 & 1) + (i >> 6 & 1) + (i >> 7 & 1) + 3;
        state_count[i] = block_count << 28 | 6;
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

