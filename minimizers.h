//
// Created by saimanalili on 25/02/2020.
//

#ifndef PH_INTEGRATION_MINIMIZERS_H
#define PH_INTEGRATION_MINIMIZERS_H

#include "common.h"

unsigned long long int getMinimizerRank(string windowSeed, int q, int windowSize);
void buildMinimizersIndexing(string& genome, string& mainName);

#endif //PH_INTEGRATION_MINIMIZERS_H
