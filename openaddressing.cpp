//
// Created by saimanalili on 16/04/2020.
//

#include "openaddressing.h"

int occupiedSpaces = 0;
int collisions = 0;

void buildOpenAddressingTables(string stringDNA, string mainName, unsigned int m, unsigned int q, unsigned long long int codeTableSize, unsigned long long int dirTableSize, unsigned long long int posTableSize) {
    auto indexing_time_start = chrono::high_resolution_clock::now();

    ofstream outfile;
    outfile.open("open_" + mainName + "_" + to_string(q) + ".txt", ios::out);

    vector<long long int> codeTable(codeTableSize);
    vector<unsigned long long int> dirTable(dirTableSize);
    vector<unsigned long long int> posTable(posTableSize);

    fill(codeTable.begin(), codeTable.end(), -1);
    fill(dirTable.begin(), dirTable.end(), 0);
    fill(posTable.begin(), posTable.end(), 0);

    vector<int> codeTableIndices;

    for (int i = 0; i < (m-q+1); i++) {
        unsigned long long int kMerIndexInGenome = extractRanking(stringDNA.substr(i, q));
        unsigned long long int hashValue = kMerIndexInGenome % (unsigned long long int)codeTableSize;
//		int hashValue = inthash_64(kMerIndexInGenome, shiftedValue - 1);

        if (codeTable.at(hashValue) != -1) {
            collisions++;
            unsigned long long int j = hashValue;
            unsigned long long int k = 1;

            while (codeTable.at(j) != -1 && codeTable.at(j) != kMerIndexInGenome) {
                j = (j + (k * k)) % (unsigned long long int)codeTableSize;
                k = k + 1;
            }

            codeTable.at(j) = kMerIndexInGenome;
            codeTableIndices.push_back(j);
            occupiedSpaces++;
        } else {
            codeTable.at(hashValue) = kMerIndexInGenome;
            codeTableIndices.push_back(hashValue);
            occupiedSpaces++;
        }
    }

    for (int i = 0; i < stringDNA.length() - q + 1; i++) {
        dirTable.at(codeTableIndices.at(i))++;
    }

    for (int i = 1; i < dirTable.size(); i++) {
        dirTable.at(i) += dirTable.at(i-1);
    }

    unsigned long long int index = 0;
    for (int i = stringDNA.length() - q; i >= 0; i--) {
        dirTable.at(codeTableIndices.at(i))--;
        index = dirTable.at(codeTableIndices.at(i));
        posTable.at(index) = i;
    }

    auto indexing_time_end = chrono::high_resolution_clock::now();
    duration<double, std::milli> indexing_duration = duration_cast<milliseconds>(indexing_time_end - indexing_time_start);
    cout << "Time taken by the indexing process is : " << (indexing_duration.count()/1000) << " sec" << endl << endl;
    cout << "Starting the writing process..." << endl << endl;

    auto writing_time_start = chrono::high_resolution_clock::now();

    outfile << "code" << endl;
    for (int i = 0; i < codeTable.size(); i++) {
        if (codeTable.at(i) != -1) {
            outfile << (unsigned long long int)codeTable.at(i) << " " << i << endl;
        } else {
            outfile << codeTable.at(i) << " " << i << endl;
        }
    }
    outfile << endl << "dir" << endl;
    for (int i = 0; i < dirTable.size(); i++) {
        outfile << dirTable.at(i) << endl;
    }
    outfile << endl << "pos" << endl;
    for (int i = 0; i < posTable.size(); i++) {
        outfile << posTable.at(i) << endl;
    }
    outfile << endl;

    auto writing_time_end = chrono::high_resolution_clock::now();
    duration<double, std::milli> writing_duration = duration_cast<milliseconds>(writing_time_end - writing_time_start);
    cout << "Time taken by the index writing process is : " << (writing_duration.count()/1000) << " sec" << endl << endl;
}

void buildOpenAddressingIndexing(string& genome, string& mainName, double loadFactor) {
    unsigned long long int codeTableSize = floor(( pow(loadFactor, -1)) * genome.size());
    unsigned long long int dirTableSize = codeTableSize + 1;
    unsigned long long int posTableSize = genome.size() - q + 1;

    cout << "Open Addressing for q = " << q  << endl << endl;
    buildOpenAddressingTables(genome, mainName, genome.length(), q, codeTableSize, dirTableSize, posTableSize);
}