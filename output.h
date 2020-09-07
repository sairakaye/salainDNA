//
// Created by saimanalili on 17/04/2020.
//

#ifndef SALAINDNA_OUTPUT_H
#define SALAINDNA_OUTPUT_H

#include "common.h"

void outputPairReads(string& mainName);
void outputRunTimeResults(string& mainName, double indexTimeTaken, double ssTimeTaken, double bmTimeTaken, double totalTimeTaken);
void outputSeedSelectorResults(string& mainName, double timeTaken);
void outputFileSeedSelectorResults(string& mainName, double timeTaken);
void outputPrealignmentResults();
void outputSAMFile();
void outputEdlibResults();

#endif //SALAINDNA_OUTPUT_H
