//
// Created by saimanalili on 17/04/2020.
//

#ifndef MULTICORE_RM_OUTPUT_H
#define MULTICORE_RM_OUTPUT_H

#include "common.h"

void outputPairReads(string& mainName);
void outputSeedSelectorResults(string& mainName, double timeTaken);
void outputFileSeedSelectorResults(string& mainName, double timeTaken);
void outputPrealignmentResults();
void outputSAMFile();

#endif //MULTICORE_RM_OUTPUT_H
