//
// Created by saimanalili on 13/04/2020.
//

#ifndef MULTICORE_RM_INDEXING_H
#define MULTICORE_RM_INDEXING_H

#include "common.h"
#include "directaddressing.h"
#include "openaddressing.h"
#include "minimizers.h"

void readIndexFile(string indexFile,
                   map<unsigned long long int, vector<unsigned long long int>>& minimizers,
                   map<long long, unsigned long long int>& codeTable,
                   vector<unsigned long long int>& dirTable,
                   vector<unsigned long long int>& posTable);

void buildIndex(string mainName, string indexFile,
                map<unsigned long long int, vector<unsigned long long int>>& minimizers,
                map<long long, unsigned long long int>& codeTable,
                vector<unsigned long long int>& dirTable,
                vector<unsigned long long int>& posTable,
                double loadFactor);

#endif //MULTICORE_RM_INDEXING_H
