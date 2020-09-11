//
// Created by saimanalili on 16/04/2020.
//

#include "openaddressing.h"

void buildOpenAddressingTables(string stringDNA, string mainName, unsigned int m, unsigned int q, unsigned long long int codeTableSize, unsigned long long int dirTableSize, unsigned long long int posTableSize)
{
    ofstream outfile;
    outfile.open( "open_" + mainName + "_" + to_string(q) + ".txt", ios::out);

    vector<long long int> initCodeTable(codeTableSize);
    vector<unsigned long long int> initDirTable(dirTableSize);
    vector<unsigned long long int> initPosTable(posTableSize);

    dirTable = vector<unsigned long long int>(initDirTable.begin(), initDirTable.end());
    initDirTable.clear();
    posTable = vector<unsigned long long int>(initPosTable.begin(), initPosTable.end());
    initPosTable.clear();

    fill(initCodeTable.begin(), initCodeTable.end(), -1);
    fill(dirTable.begin(), dirTable.end(), 0);
    fill(posTable.begin(), posTable.end(), 0);

    vector<unsigned long long int> codeTableIndices(m - q + 1);
    int i;

    omp_lock_t myLock;
    omp_init_lock(&myLock);

    #pragma omp parallel for
    for (i = 0; i < (m - q + 1); i++) {
        unsigned long long int kMerIndexInGenome = extractRanking(stringDNA.substr(i, q));
        unsigned long long int hashValue = kMerIndexInGenome % (unsigned long long int)codeTableSize;

        omp_set_lock(&myLock);
        if (initCodeTable.at(hashValue) != -1) {
            unsigned long long int j = hashValue;
            unsigned long long int k = 1;

            while (initCodeTable.at(j) != -1 && initCodeTable.at(j) != kMerIndexInGenome) {
                j = (j + (k * k)) % (unsigned long long int) codeTableSize;
                k = k + 1;
            }
            initCodeTable.at(j) = kMerIndexInGenome;
            dirTable.at(j)++;

            omp_unset_lock(&myLock);

            codeTableIndices.at(i) = j;
        } else {
            initCodeTable.at(hashValue) = kMerIndexInGenome;
            dirTable.at(hashValue)++;

            omp_unset_lock(&myLock);

            codeTableIndices.at(i) = hashValue;
        }
    }

    omp_destroy_lock(&myLock);

    for (i = 1; i < dirTable.size(); i++) {
        dirTable.at(i) += dirTable.at(i-1);
    }

    #pragma omp parallel for
    for (i = stringDNA.length() - q; i >= 0; i--) {
        unsigned long long int index = 0;

        #pragma omp critical
        {
            dirTable.at(codeTableIndices.at(i))--;
            index = dirTable.at(codeTableIndices.at(i));
        }
        posTable.at(index) = i;
    }

    outfile << "code" << endl;
    i = 0;
    auto it = initCodeTable.begin();
    while (it != initCodeTable.end())
    {
        auto currIt = it++;

        if (*currIt != -1) {
            outfile << (unsigned long long int)*currIt << " " << i << endl;
            //codeTable[i] = (unsigned long long int)*currIt;

            codeTable.insert(pair<long long, unsigned long long int>(*currIt, i));
        } else {
            outfile << *currIt << " " << i << endl;
            //codeTable[i] = *currIt;
        }



//        initCodeTable.erase(currIt);
        i++;
    }

    /*
    for (int i = 0; i < initCodeTable.size(); i++) {
        if (initCodeTable.at(i) != -1) {
            outfile << (unsigned long long int)initCodeTable.at(i) << " " << i << endl;
        } else {
            outfile << initCodeTable.at(i) << " " << i << endl;
        }

        //codeTable[i] = (unsigned long long int)initCodeTable.at(i);
    }
    */

    outfile << endl << "dir" << endl;
    for (int i = 0; i < dirTable.size(); i++) {
        outfile << dirTable.at(i) << endl;
    }
    outfile << endl << "pos" << endl;
    for (int i = 0; i < posTable.size(); i++) {
        outfile << posTable.at(i) << endl;
    }
    outfile << endl;

    outfile.close();
}

void buildOpenAddressingIndexing(string& genome, string& mainName, double loadFactor) {
    unsigned long long int codeTableSize = floor(( pow(loadFactor, -1)) * genome.size());
    unsigned long long int dirTableSize = codeTableSize + 1;
    unsigned long long int posTableSize = genome.size() - q + 1;

    cout << "Open Addressing for q = " << q  << endl << endl;
    buildOpenAddressingTables(genome, mainName, genome.length(), q, codeTableSize, dirTableSize, posTableSize);
}