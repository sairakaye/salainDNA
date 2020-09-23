/**
 * @file directaddressing.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the implementation of the Direct Addressing mode of Hash-based Indexer.
 */

#include "directaddressing.h"

/**
 * It builds the direct addressing index.
 *
 * @param genome - The reference genome data.
 * @param genomeFileName - The filename of the reference genome.
 * @param n - The length of the reference genome.
 * @param q - The value of the q-gram.
 * @param dirTableSize - The length or size of the directory table.
 * @param posTableSize - The length or size of the position table.
 */
void buildDirectAddressingTables(string genome, string genomeFileName, unsigned int n, unsigned int q, unsigned long long int dirTableSize, unsigned long long int posTableSize)
{
    ofstream outfile;
    outfile.open("dir_" + genomeFileName + "_" + to_string(q) + ".txt", ios::out);

    vector<unsigned long long int> initDirTable(dirTableSize);
    vector<unsigned long long int> initPosTable(posTableSize);

    dirTable = vector<unsigned long long int>(initDirTable.begin(), initDirTable.end());
    initDirTable.clear();
    posTable = vector<unsigned long long int>(initPosTable.begin(), initPosTable.end());
    initPosTable.clear();

    fill(dirTable.begin(), dirTable.end(), 0);
    fill(posTable.begin(), posTable.end(), 0);

    int i;

    #pragma omp parallel for
    for (i = 0; i < n - q + 1; i++) {
        unsigned long long int rank = extractRanking(genome.substr(i, q));

        #pragma omp atomic
        dirTable.at(rank)++;
    }

    for (i = 1; i < dirTable.size(); i++) {
        dirTable.at(i) += dirTable.at(i-1);
    }

    #pragma omp parallel for
    for (i = n - q; i >= 0; i--) {
        unsigned long long int rank = extractRanking(genome.substr(i, q));
        unsigned long long int index = 0;

        #pragma omp critical
        {
            dirTable.at(rank)--;
            index = dirTable.at(rank);
        }
        posTable.at(index) = i;
    }

    outfile << "dir" << "\n";
    for (int i = 0; i < dirTable.size(); i++) {
        outfile << dirTable.at(i) << "\n";
    }
    outfile << "\n" << "pos" << "\n";
    for (int i = 0; i < posTable.size(); i++) {
        outfile << posTable.at(i) << "\n";
    }

    outfile.close();
}

/**
 * It initializes the table sizes and calls the building of the direct addressing index.
 *
 * @param genome - The reference genome data.
 * @param genomeFileName - The filename of the reference genome.
 */
void buildDirectAddressingIndexingFile(string& genome, string& genomeFileName) {
    unsigned long long int dirTableSize = pow(4, q) + 1;
    unsigned long long int posTableSize = genome.size() - q + 1;

    cout << "Direct Addressing for q = " << q << endl << endl;
    buildDirectAddressingTables(genome, genomeFileName, genome.length(), q, dirTableSize, posTableSize);
}