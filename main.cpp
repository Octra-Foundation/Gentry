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

template <class Template> void alloc(Template*& get_next_bit, int bit_count) {
    get_next_bit = (Template*)calloc(bit_count, sizeof(Template));
    if (!get_next_bit) {
        static_cast<void>(fprintf(stderr, "out of memory\n")), exit(1);
    }
}

