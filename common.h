//
// Created by saimanalili on 25/02/2020.
//

#ifndef PH_INTEGRATION_COMMON_H
#define PH_INTEGRATION_COMMON_H

#include <iostream>
#include <fstream>
#include <map>
#include <iterator>
#include <string>
#include <utility>
#include <algorithm>
#include <vector>
#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

extern string refGenome;
extern vector<string> reads;

extern string mode;
extern string searchMode;
extern unsigned int q;
extern unsigned int w;
extern unsigned int m;
extern unsigned int e;

extern map<unsigned long long int, vector<unsigned long long int>> minimizers;
extern map<long long, unsigned long long int> codeTable;
extern vector<unsigned long long int> dirTable;
extern vector<unsigned long long int> posTable;

extern map<string, vector<unsigned long long int>> forwardReadsMap;
extern map<string, vector<unsigned long long int>> reverseReadsMap;

typedef struct {
    unsigned int numSeeds;
    unsigned int numReads;
    unsigned int numAcceptedSeeds;
    unsigned int numAcceptedReads;
} Counters;

extern Counters *counter;

extern string temp_comp;
extern ofstream outputPossibleReadsFile;
extern ofstream outputLocationsFile;
extern ofstream infoFile;

string readGenomeFile(string filename);
vector<string> readReadsFile(string filename);
void getDirectAddressing(string filename, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable);
void getOpenAddressing(string filename, map<long long, unsigned long long int>& codeTable, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable);
map<unsigned long long int, vector<unsigned long long int>> getMinimizers(string filename);
void results(map<string, vector<unsigned long long int>>& forwardReadsMap, map<string, vector<unsigned long long int>>& reverseReadsMap);

#endif //PH_INTEGRATION_COMMON_H
