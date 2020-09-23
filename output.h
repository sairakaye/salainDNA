/**
 * @file output.h
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief Header file for output of the results, containing function declarations to be used by the program.
 */

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

#endif
