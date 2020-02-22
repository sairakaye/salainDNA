//
// Created by saimanalili on 12/02/2020.
//

#ifndef PIGEONHOLE_PHREADS_H
#define PIGEONHOLE_PHREADS_H

#include "pigeonhole.h"
#include "minimizers.h"

void partitioningReadsToSeeds(string filename, string mode, int windowLength, int q);
void partitioningReadsToSeedsExit(string filename, string mode, int windowLength, int q);

#endif //PIGEONHOLE_PHREADS_H
