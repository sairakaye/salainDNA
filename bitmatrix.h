//
// Created by saimanalili on 29/03/2020.
//

#ifndef PIGEONHOLE_BITMATRIX_BITMATRIX_H
#define PIGEONHOLE_BITMATRIX_BITMATRIX_H

#include "common.h"

#include <zconf.h>

extern int TruePos;
extern int TrueNeg;
extern int FalsePos;
extern int FalseNeg;

extern int alignmentNeeded;
extern int notNeeded;

//void multiThreadedMain(map<string, vector<unsigned long long int>>& forwardReadsMap, map<string, vector<unsigned long long int>>& reverseReadsMap);
void multiThreadedMain();

#endif //PIGEONHOLE_BITMATRIX_BITMATRIX_H
