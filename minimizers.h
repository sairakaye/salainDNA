//
// Created by saimanalili on 25/02/2020.
//

#ifndef SALAINDNA_MINIMIZERS_H
#define SALAINDNA_MINIMIZERS_H

#include "common.h"

unsigned long long int getMinimizerRank(string windowSeed, int q, int windowSize);
unsigned long long int getMinimizerRankWithoutWindow(string windowSeed, int q);
void buildMinimizersIndexing(string& genome, string& mainName);

#endif //PH_INTEGRATION_MINIMIZERS_H
