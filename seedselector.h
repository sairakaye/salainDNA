/**
 * @file seedselector.h
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief Header file for Seed Selector, containing variables and function declarations to be used by the program.
 */

#ifndef SALAINDNA_PIGEONHOLE_H
#define SALAINDNA_PIGEONHOLE_H

#include "common.h"
#include "searchall.h"
#include "searchexit.h"
#include "searching.h"

extern int j; // The number of partitions in each read.
extern int allowableE; // The value of allowable error threshold in each seed.

void seedSelectorSearchAllProcess();
void seedSelectorExitProcess();

#endif
