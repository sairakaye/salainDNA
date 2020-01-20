//
// Created by saimanalili on 19/01/2020.
//

#ifndef PIGEONHOLE_ADDRESSING_H
#define PIGEONHOLE_ADDRESSING_H

#include "common.h"
#include "minimizers.h"

void buildTablesDirect(string stringDNA, int q, double dirTableSize, double posTableSize, vector<int> dirTable, vector<int> posTable);
void buildTablesOpen(string stringDNA, int q, int shiftedValue, double codeTableSize, double dirTableSize, double posTableSize, vector<int> codeTable, vector<int> dirTable, vector<int> posTable);
#endif //PIGEONHOLE_ADDRESSING_H
