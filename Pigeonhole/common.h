//
// Created by saimanalili on 29/12/2019.
//
#ifndef MP_PIGEONHOLE_COMMON_H
#define MP_PIGEONHOLE_COMMON_H

#include <iostream>
#include <fstream>
#include <map>
#include <iterator>
#include <string>
#include <utility>
#include <algorithm>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

extern string refGenome;
extern vector<string> seeds;
extern map<unsigned long int, vector<unsigned long int>> minimizers;

extern map<long long, unsigned long int> codeTable;
extern vector<unsigned long int> dirTable;
extern vector<unsigned long int> posTable;

extern vector<unsigned long int> forwardFound;
extern vector<unsigned long int> reverseFound;

typedef struct {
    string kmer;
    int rank;
} Code;

//extern Indexing* indexing;
//extern DirectIndexing *directIndexing;

//typedef struct {
//    char *genome;
//    int length;
//} RefGenome;
//
//typedef struct {
//    char **seeds;
//    int size;
//} ReadList;

string getFilename(string filename);
string readGenome(string filename);
vector<string> readReads(string filename);
map<unsigned long int, vector<unsigned long int>> getMinimizersFromFile(string filename);
//RefGenome *readGenome(std::string filename);
//ReadList *readReads(std::string filename);
void generateQGrams(string prefix, vector<string>& qGrams, int k);
void getDirectAddressing(string filename, vector<unsigned long int>& dirTable, vector<unsigned long int>& posTable);
void getOpenAddressing(string filename, map<long long, unsigned long int>& codeTable, vector<unsigned long int>& dirTable, vector<unsigned long int>& posTable);

#endif //MP_PIGEONHOLE_COMMON_H
