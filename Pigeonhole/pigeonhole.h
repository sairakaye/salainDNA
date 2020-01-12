//
// Created by saimanalili on 30/12/2019.
//

#ifndef MP_PIGEONHOLE_PIGEONHOLE_H
#define MP_PIGEONHOLE_PIGEONHOLE_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>

#include "indexing.h"

typedef struct {
    char *read;
    char *readFromGenome;
    //std::string read;
    //std::string readFromGenome;
} AlignedReads;

void filterReadsWithMinimizers(std::string filename, ReadList *readList, std::map<int, std::vector<int>> minimizers, RefGenome *refGenome, int e, int segmentLength);
void filterReads(std::string filename, ReadList *readList, HashIndexing *hashIndexing, RefGenome *refGenome, int e, int segmentLength);
void parallelizeFilterReads(std::string filename, ReadList *readList, HashIndexing *hashIndexing, RefGenome *refGenome, int e, int segmentLength) ;
#endif //MP_PIGEONHOLE_PIGEONHOLE_H