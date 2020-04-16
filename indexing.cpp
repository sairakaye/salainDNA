//
// Created by saimanalili on 13/04/2020.
//

#include "indexing.h"

void readIndexFile(string mode, string indexFile,
                   map<unsigned long long int, vector<unsigned long long int>>& minimizers,
                   map<long long, unsigned long long int>& codeTable,
                   vector<unsigned long long int>& dirTable,
                   vector<unsigned long long int>& posTable) {
    cout << "Reading the indexing... " << endl << indexFile << endl << endl;

    if (mode.compare("min") == 0) {
        minimizers = getMinimizers(indexFile);
    } else if (mode.compare("dir") == 0) {
        getDirectAddressing(indexFile, dirTable, posTable);
    } else if (mode.compare("open") == 0) {
        getOpenAddressing(indexFile, codeTable, dirTable, posTable);
    } else {
        cout << "Mode not valid...";
        exit(EXIT_FAILURE);
    }
}

void buildIndex(string mode, string mainName, string indexFile,
                map<unsigned long long int, vector<unsigned long long int>>& minimizers,
                map<long long, unsigned long long int>& codeTable,
                vector<unsigned long long int>& dirTable,
                vector<unsigned long long int>& posTable) {
    cout << "Starting the indexing... " << endl << indexFile << endl << endl;

    if (mode.compare("min") == 0) {
        buildMinimizersIndexing(refGenome, mainName);
    } else if (mode.compare("dir") == 0) {
        buildDirectAddressingIndexing(refGenome, mainName);
    } else if (mode.compare("open") == 0) {
        buildOpenAddressingIndexing(refGenome, mainName);
    } else {
        cout << "Mode not valid...";
        exit(EXIT_FAILURE);
    }
}

