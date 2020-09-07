//
// Created by saimanalili on 29/03/2020.
//

#ifndef SALAINDNA_BITMATRIX_H
#define SALAINDNA_BITMATRIX_H

#include "common.h"

//#include <zconf.h>


extern int truePos;
extern int trueNeg;

/*
extern int FalsePos;
extern int FalseNeg;
*/

extern unsigned int alignmentNeeded;
extern unsigned int notNeeded;

//void bitMatrixFilterProcess(map<string, vector<unsigned long long int>>& forwardReadsMap, map<string, vector<unsigned long long int>>& reverseReadsMap);
void bitMatrixFilterProcess();
void verifyWithEdlib();
void preCheckWithEdlib();

#endif //SALAINDNA_BITMATRIX_H
