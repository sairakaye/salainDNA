//
// Created by saimanalili on 25/02/2020.
//

#ifndef PH_INTEGRATION_SEARCHING_H
#define PH_INTEGRATION_SEARCHING_H

#include "common.h"
#include "minimizers.h"

void exactSearchingPosition(string seed, string mode, int q, int k, vector<unsigned long long int>& foundLocations);
void approximateSearchingPosition(string seed, string mode, int q, int k, vector<unsigned long long int>& foundLocations, int *minEditFound);

#endif //PH_INTEGRATION_SEARCHING_H
