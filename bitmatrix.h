//
// Created by saimanalili on 29/03/2020.
//

#ifndef MULTICORE_RM_BITMATRIX_H
#define MULTICORE_RM_BITMATRIX_H

#include "common.h"

//#include <zconf.h>

/*
extern int TruePos;
extern int TrueNeg;
extern int FalsePos;
extern int FalseNeg;
*/

extern unsigned int alignmentNeeded;
extern unsigned int notNeeded;

//void multiThreadedMain(map<string, vector<unsigned long long int>>& forwardReadsMap, map<string, vector<unsigned long long int>>& reverseReadsMap);
void multiThreadedMain();

#endif //MULTICORE_RM_BITMATRIX_H
