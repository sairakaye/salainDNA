//
// Created by saimanalili on 25/02/2020.
//

#ifndef PH_INTEGRATION_SEARCHING_H
#define PH_INTEGRATION_SEARCHING_H

#include "common.h"
#include "minimizer.h"
#include "edlib.h"

void searchingPosition(string seed, string read, string mode, int q, int k, bool isForwardStrand, bool isCheckApproximate, vector<unsigned long long int>& foundLocations);

#endif //PH_INTEGRATION_SEARCHING_H
