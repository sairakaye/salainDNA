//
// Created by saimanalili on 17/04/2020.
//

#ifndef SALAINDNA_OUTPUT_H
#define SALAINDNA_OUTPUT_H

#include "common.h"

void outputPairReads(string& genomeFileName);
void outputRunTimeResults(string& genomeFileName, double indexRunTime, double ssRunTime, double bmRunTime,
                          double verificationRunTime, double totalRunTime);
void outputSeedSelectorResults(string& genomeFileName, double timeTaken);
void outputPrealignmentResults();
void outputSAMFile();
void outputEdlibResults();

#endif //SALAINDNA_OUTPUT_H
