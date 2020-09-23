/**
 * @file bitmatrix.h
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief Header file for Bit Matrix, containing variables and function declarations to be used by the program.
 */

#ifndef SALAINDNA_BITMATRIX_H
#define SALAINDNA_BITMATRIX_H

#include "common.h"

extern int truePos; // The number  of True Positive locations filtered.
extern int trueNeg; // The number of True Negative locations filtered.

extern unsigned int alignmentNeeded; // The number of accepted filtered locations.
extern unsigned int notNeeded; // The number of rejected filtered locations.

void bitMatrixFilterProcess();

#endif
