//
// Created by saimanalili on 25/02/2020.
//

#ifndef SALAINDNA_COMMON_H
#define SALAINDNA_COMMON_H

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

typedef struct {
    string genomeName;
    string genomeData;
} Genome;

typedef struct {
    string readName;
    string readData;
    vector<unsigned long long int> forwardLocations;
    vector<unsigned long long int> reverseLocations;
} Read;


extern vector<pair<string, int>> alphabetRef;
extern Genome refGenome;
extern vector<Read> reads;

extern string mode;
extern string searchMode;
extern bool isReverseAccepted;
extern unsigned int q;
extern unsigned int w;
extern unsigned int m;
extern unsigned int e;
extern double loadFactor;

extern string SAMFileName;

extern map<unsigned long long int, vector<unsigned long long int>> minimizers;
extern map<long long, unsigned long long int> codeTable;
extern vector<unsigned long long int> dirTable;
extern vector<unsigned long long int> posTable;

extern unsigned int numSeeds;
extern unsigned int numReads;
extern unsigned int numAcceptedSeeds;
extern unsigned int numAcceptedReads;
extern unsigned int numPossibleReadLocations;
extern unsigned int numFilteredReadLocations;
extern unsigned int numVerifiedReadLocations;

extern double indexRunTime;
extern double ssRunTime;
extern double bmRunTime;
extern double verificationRunTime;

unsigned long long int extractRanking(string kMer);
uint64_t inthash_64(uint64_t key, uint64_t mask);
Genome readGenomeFile(string filename);
void readReadsFile(string filename);
string reverseComplement(string read);
void getDirectAddressing(string filename, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable);
void getOpenAddressing(string filename, map<long long, unsigned long long int>& codeTable, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable);
map<unsigned long long int, vector<unsigned long long int>> getMinimizers(string filename);

#endif //SALAINDNA_COMMON_H
