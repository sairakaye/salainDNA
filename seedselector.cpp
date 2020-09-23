/**
 * @file seedselector.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the Seed Selector implementation.
 */


#include "seedselector.h"

int j; // The number of partitions in each read.
int allowableE; // The value of allowable error threshold in each seed.

/**
* Starts the seed selector process using the search all mode.
*
*/
void seedSelectorSearchAllProcess() {
    j = ceil(m / (double) q);
    allowableE = floor(e / (double) j);

    if (allowableE == 0) {
        exactSearchingForAll();
    } else {
        approximateSearchingForAll();
    }
}

/**
* Starts the seed selector process using the exit mode.
*
*/
void seedSelectorExitProcess() {
    j = ceil(m / (double) q);
    allowableE = floor(e / (double) j);

    if (allowableE == 0) {
        exactSearchingForExit();
    } else {
        approximateSearchingForExit();
    }
}