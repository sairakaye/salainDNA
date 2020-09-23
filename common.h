/**
 * @file common.h
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief Acted as the main header file, containing all global variables and functions needed in the program.
 */

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

/**
 * Declaration of a struct Genome. This contains the data of the genome
 * as well as the name of the genome if specified in the file.
 */
typedef struct {
    string genomeName;
    string genomeData;
} Genome;

/**
 * Declaration of a struct Read. This contains the data of the read,
 * the name of the read if specified in the file, and the locations
 * both from when it is in forward strand and backward strand.
 */
typedef struct {
    string readName;
    string readData;
    vector<unsigned long long int> forwardLocations;
    vector<unsigned long long int> reverseLocations;
} Read;

extern vector<pair<string, int>> alphabetRef; // It contains the equivalent numerical value of each DNA alphabet code.
extern Genome refGenome; // Data of the reference genome.
extern vector<Read> reads; // A list of reads to be mapped by SalainDNA.

extern string mode; // It indicates what mode of hash-based indexer will be used.
extern string searchMode; // It indicates what searching mode of the seed selector will be used.
extern bool isReverseAccepted; // It indicates if the reverse complement of each read will also be mapped or not.
extern unsigned int q; // The value of the q-gram.
extern unsigned int w; // The window value to be used in minimizers hash-based indexer.
extern unsigned int m; // The length of each read.
extern unsigned int e; // The value of error threshold.
extern double loadFactor; // The value of load factor used by the code table from the open addressing hash-based indexer.

extern string SAMFileName; // The file name of the output SAM file from SalainDNA.

extern map<unsigned long long int, vector<unsigned long long int>> minimizers; // The data structure for minimizer-based indexing.
extern map<long long, unsigned long long int> codeTable; // It is where q-gram ranks are hashed into.
extern vector<unsigned long long int> dirTable; // It contains the starting location of q-grams in the position table.
extern vector<unsigned long long int> posTable; // It is a list containing the positions of q-grams in the reference genome.

extern unsigned int numSeeds; // Overall number of seeds checked by SalainDNA.
extern unsigned int numReads; // Overall number of reads checked by SalainDNA.
extern unsigned int numAcceptedSeeds; // Overall number of accepted seeds by SalainDNA.
extern unsigned int numAcceptedReads; // Overall number of accepted reads by SalainDNA.
extern unsigned int numPossibleReadLocations; // Number of possible read locations found from Seed Selector.
extern unsigned int numFilteredReadLocations; // Number of filtered read locations from Bit Matrix.
extern unsigned int numVerifiedReadLocations; // Number of verified read locations from Edlib.

extern double indexRunTime; // The runtime of the Hash-based Indexer.
extern double ssRunTime; // The runtime of Seed Selector.
extern double bmRunTime; // The runtime of the Bit Matrix.
extern double verificationRunTime; // The runtime of Edlib.

unsigned long long int extractRanking(string kMer);
uint64_t inthash_64(uint64_t key, uint64_t mask);
Genome readGenomeFile(string filename);
vector<Read> readReadsFile(string filename);
string reverseComplement(string read);
void getDirectAddressing(string filename, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable);
void getOpenAddressing(string filename, map<long long, unsigned long long int>& codeTable, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable);
map<unsigned long long int, vector<unsigned long long int>> getMinimizers(string filename);

#endif
