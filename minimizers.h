/**
 * @file minimizers.h
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief Header file for the Minimizers, containing function declarations to be used by the program.
 */

#ifndef SALAINDNA_MINIMIZERS_H
#define SALAINDNA_MINIMIZERS_H

#include "common.h"

unsigned long long int getMinimizerRank(string windowSeed, int q, int windowSize);
unsigned long long int getMinimizerRankWithoutWindow(string windowSeed, int q);
void buildMinimizersIndexingFile(string& genome, string& genomeFileName);

#endif
