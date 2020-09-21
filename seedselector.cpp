//
// Created by saimanalili on 25/02/2020.
//

#include "seedselector.h"

int j;
int allowableE;

void seedSelectorSearchAllProcess() {
    j = ceil(m / (double) q);
    allowableE = floor(e / (double) j);

    if (allowableE == 0) {
        exactSearchingForAll();
    } else {
        approximateSearchingForAll();
    }
}

void seedSelectorExitProcess() {
    j = ceil(m / (double) q);
    allowableE = floor(e / (double) j);

    if (allowableE == 0) {
        exactSearchingForExit();
    } else {
        approximateSearchingForExit();
    }
}