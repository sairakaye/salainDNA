//
// Created by saimanalili on 25/02/2020.
//

#ifndef MULTICORE_RM_COMMON_H
#define MULTICORE_RM_COMMON_H

#include <iostream>
#include <fstream>
#include <map>
#include <iterator>
#include <string>
#include <utility>
#include <algorithm>
#include <vector>
#include <bits/stdc++.h>
#include <chrono>
#include <omp.h>

#include "edlib.h"

using namespace std;
using namespace std::chrono;

extern string refGenome;
//extern vector<string> reads;
extern map<string, string> reads;

extern string mode;
extern string searchMode;
extern unsigned int q;
extern unsigned int w;
extern unsigned int m;
extern unsigned int e;
extern double loadFactor;

extern vector<pair<string, int> > alphabetRef;

extern map<unsigned long long int, vector<unsigned long long int>> minimizers;
extern map<long long, unsigned long long int> codeTable;
extern vector<unsigned long long int> dirTable;
extern vector<unsigned long long int> posTable;

extern map<string, vector<string>> readsLabelMap;
extern map<string, vector<unsigned long long int>> forwardReadsMap;
extern map<string, vector<unsigned long long int>> reverseReadsMap;
extern map<string, vector<unsigned long long int>> filteredReadsMap;

extern unsigned int numSeeds;
extern unsigned int numReads;
extern unsigned int numAcceptedSeeds;
extern unsigned int numAcceptedReads;
extern unsigned int numLocationsForward;
extern unsigned int numLocationsReverse;

unsigned long long int extractRanking(string kMer);
uint64_t inthash_64(uint64_t key, uint64_t mask);
string readGenomeFile(string filename);
map<string, string> readReadsFile(string filename);
string reverseComplement(string read);
void getDirectAddressing(string filename, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable);
void getOpenAddressing(string filename, map<long long, unsigned long long int>& codeTable, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable);
map<unsigned long long int, vector<unsigned long long int>> getMinimizers(string filename);
void removingDuplicateLocationsInEachRead();

#endif //MULTICORE_RM_COMMON_H
