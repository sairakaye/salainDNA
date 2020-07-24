//
// Created by saimanalili on 16/04/2020.
//

#include "directaddressing.h"

void buildDirectAddressingTables(string stringDNA, string mainName, unsigned int m, unsigned int q, unsigned long long int dirTableSize, unsigned long long int posTableSize)
{
    ofstream outfile;
    outfile.open("dir_" + mainName + "_" + to_string(q) + ".txt", ios::out);

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
    for (i = 0; i < (m - q + 1); i++) {
        unsigned long long int tempIndex = extractRanking(stringDNA.substr(i, q));

        #pragma omp atomic
        dirTable.at(tempIndex)++;
    }

    for (i = 1; i < dirTable.size(); i++) {
        dirTable.at(i) += dirTable.at(i-1);
    }

    #pragma omp parallel for
    for (i = m - q; i >= 0; i--) {
        unsigned long long int tempIndex = extractRanking(stringDNA.substr(i, q));
        unsigned long long int index = 0;

        #pragma omp critical
        {
            dirTable.at(tempIndex)--;
            index = dirTable.at(tempIndex);
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

void buildDirectAddressingIndexing(string& genome, string& mainName) {
    unsigned long long int dirTableSize = pow(4, q) + 1;
    unsigned long long int posTableSize = genome.size() - q + 1;

    cout << "Direct Addressing for q = " << q << endl << endl;
    buildDirectAddressingTables(genome, mainName, genome.length(), q, dirTableSize, posTableSize);
}