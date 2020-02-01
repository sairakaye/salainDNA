//
// Created by saimanalili on 08/01/2020.
//

#ifndef MP_PIGEONHOLE_MINIMIZERS_H
#define MP_PIGEONHOLE_MINIMIZERS_H

#include "common.h"

extern vector<pair<string, int> > alphabetRef;

unsigned long int extractRanking(string kmer);
uint64_t inthash_64(uint64_t key, uint64_t mask);
unsigned int getMinimizerRank(string windowSeed, int q, int windowSize);
//multimap<vector<int>, int> generateMinimizers(string stringDNA, int q, int w, int m);

#endif //MP_PIGEONHOLE_MINIMIZERS_H
