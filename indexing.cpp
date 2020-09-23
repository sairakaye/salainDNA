/**
 * @file indexing.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the implementation of the Hash-based indexer.
 */

#include "indexing.h"

/**
 * It reads the index file indicated in the program.
 *
 * @param indexFile - The filename of the index file.
 * @param minimizers - It is a map containing locations of each minimizer rank.
 * @param codeTable - It is where q-gram ranks are hashed into.
 * @param dirTable - It contains the starting location of q-grams in the position table.
 * @param posTable - It is a list containing the positions of q-grams in the reference genome.
 */
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

/**
 * It builds the hash-based index specified in the read mapper.
 *
 * @param genomeFileName - The filename of the reference genome.
 * @param indexFile - The filename of the index file.
 * @param minimizers - It is a map containing locations of each minimizer rank.
 * @param codeTable - It is where q-gram ranks are hashed into.
 * @param dirTable - It contains the starting location of q-grams in the position table.
 * @param posTable - It is a list containing the positions of q-grams in the reference genome.
 * @param loadFactor - The value of load factor used by the code table from the open addressing hash-based indexer.
 */
void buildIndex(string genomeFileName, string indexFile,
                map<unsigned long long int, vector<unsigned long long int>>& minimizers,
                map<long long, unsigned long long int>& codeTable,
                vector<unsigned long long int>& dirTable,
                vector<unsigned long long int>& posTable,
                double loadFactor) {
    cout << "Generating the index file... " << endl << indexFile << endl << endl;

    if (mode.compare("min") == 0) {
        buildMinimizersIndexingFile(refGenome.genomeData, genomeFileName);
    } else if (mode.compare("dir") == 0) {
        buildDirectAddressingIndexingFile(refGenome.genomeData, genomeFileName);
    } else if (mode.compare("open") == 0) {
        buildOpenAddressingIndexingFile(refGenome.genomeData, genomeFileName, loadFactor);
    } else {
        cout << "Mode not valid...";
        exit(EXIT_FAILURE);
    }
}

