//
// Created by saimanalili on 13/04/2020.
//

#include "indexing.h"

void readIndexFile(string indexFile,
                   map<unsigned long long int, vector<unsigned long long int>>& minimizers,
                   map<long long, unsigned long long int>& codeTable,
                   vector<unsigned long long int>& dirTable,
                   vector<unsigned long long int>& posTable) {
    cout << "Reading the index file... " << endl << indexFile << endl << endl;

    if (mode.compare("min") == 0) {
        minimizers = getMinimizers(indexFile);
    } else if (mode.compare("dir") == 0) {
        getDirectAddressing(indexFile, dirTable, posTable);

        if (posTable.size() < refGenome.genomeData.length() - q + 1) {
            cout << "Mismatch size of position table." << endl;
            exit(EXIT_FAILURE);
        }
    } else if (mode.compare("open") == 0) {
        getOpenAddressing(indexFile, codeTable, dirTable, posTable);

        if (posTable.size() < refGenome.genomeData.length() - q + 1) {
            cout << "Mismatch size of position table." << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        cout << "Mode not valid..." << endl;
        exit(EXIT_FAILURE);
    }
}

void buildIndex(string mainName, string indexFile,
                map<unsigned long long int, vector<unsigned long long int>>& minimizers,
                map<long long, unsigned long long int>& codeTable,
                vector<unsigned long long int>& dirTable,
                vector<unsigned long long int>& posTable,
                double loadFactor) {
    cout << "Generating the index file... " << endl << indexFile << endl << endl;

    if (mode.compare("min") == 0) {
        buildMinimizersIndexing(refGenome.genomeData, mainName);
    } else if (mode.compare("dir") == 0) {
        buildDirectAddressingIndexing(refGenome.genomeData, mainName);
    } else if (mode.compare("open") == 0) {
        buildOpenAddressingIndexing(refGenome.genomeData, mainName, loadFactor);
    } else {
        cout << "Mode not valid...";
        exit(EXIT_FAILURE);
    }
}

