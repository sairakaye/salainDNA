/**
 * @file searching.h
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief Header file for location searching of Seed Selector, containing function declarations to be used by the program.
 */

#ifndef SALAINDNA_SEARCHING_H
#define SALAINDNA_SEARCHING_H

#include "common.h"
#include "minimizers.h"

vector<unsigned long long int> exactSearchingPosition(string seed, string mode, int adjustmentValue);
vector<unsigned long long int> approximateSearchingPosition(string seed, string mode, int adjustmentValue);

#endif //SALAINDNA_SEARCHING_H
