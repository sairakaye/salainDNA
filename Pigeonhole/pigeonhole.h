//
// Created by saimanalili on 30/12/2019.
//

#ifndef MP_PIGEONHOLE_PIGEONHOLE_H
#define MP_PIGEONHOLE_PIGEONHOLE_H

#include "common.h"

typedef struct {
   string read;
   string readFromGenome;
} AlignedReads;


void filterReads(string filename, int q, int m, int j, int k);
void parallelizeFilterReads(string filename, int q, int m, int j, int k);

#endif //MP_PIGEONHOLE_PIGEONHOLE_H