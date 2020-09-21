//
// Created by saimanalili on 25/02/2020.
//

#ifndef SALAINDNA_SEARCHING_H
#define SALAINDNA_SEARCHING_H

#include "common.h"
#include "minimizers.h"

vector<unsigned long long int> exactSearchingPosition(string seed, string mode, int adjustmentValue);
vector<unsigned long long int> approximateSearchingPosition(string seed, string mode, int adjustmentValue);

#endif //SALAINDNA_SEARCHING_H
