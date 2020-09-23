/**
 * @file indexing.h
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief Header file for the Hash-based Indexer, containing function declarations to be used by the program.
 */

#ifndef SALAINDNA_INDEXING_H
#define SALAINDNA_INDEXING_H

#include "common.h"
#include "directaddressing.h"
#include "openaddressing.h"
#include "minimizers.h"

void readIndexFile(string indexFile,
                   map<unsigned long long int, vector<unsigned long long int>>& minimizers,
                   map<long long, unsigned long long int>& codeTable,
                   vector<unsigned long long int>& dirTable,
                   vector<unsigned long long int>& posTable);

void buildIndex(string genomeFileName, string indexFile,
                map<unsigned long long int, vector<unsigned long long int>>& minimizers,
                map<long long, unsigned long long int>& codeTable,
                vector<unsigned long long int>& dirTable,
                vector<unsigned long long int>& posTable,
                double loadFactor);

#endif
