#include "common.h"
#include "pigeonhole.h"
#include "addressing.h"

string refGenome;
vector<string> seeds;
map<unsigned long long, vector<unsigned long long>> minimizers;
vector<string> qgrams;
vector<int> codeTable;
vector<int> dirTable;
vector<int> posTable;

int main(int argc, char *argv[]) {
    cout << "Reading the reference genome..." << endl;

    string filepath = "/home/saimanalili/multicore-rm/Pigeonhole/genome/sim_1m.fa";

    filepath.erase(remove(filepath.begin(), filepath.end(), '\"'), filepath.end());
    filepath.erase(remove(filepath.begin(), filepath.end(), '\''), filepath.end());

    refGenome = readGenome(filepath);

    //refGenome = "TTATCTCTTA";

    cout << "Done!" << endl;

    /**
     * m - length of the read.
     * q - length of the seed (k-mer/q-gram)
     * k - errors
     * j - number of partitions in the read.
     */

    int m = 100;
    int q = 14;
    int w = 14;
    int k = 2;
    int j = m / q;
    int windowLength = w + q - 1;
    //float loadFactor = 0.8;

    /**
     * For Minimizers.
     */
    cout << "Processing minimizers..." << endl;
    string fileMinimizers = "/home/saimanalili/multicore-rm/Pigeonhole/minimizers/min_sim_1m_14_inthash.txt";
    minimizers = getMinimizersFromFile(fileMinimizers);
    cout << "Done!" << endl;


    /**
     * For Direct Addressing
     */
//    generateQGrams("", qgrams, q);
//    getDirectAddressing("mama", qgrams, dirTable, posTable);


    cout << "Reading the seeds..." << endl;
    string readsFilename = "/home/saimanalili/multicore-rm/Pigeonhole/reads/random_27bp_1000r.fa";
    seeds = readReads(readsFilename);
    cout << "Done!" << endl;

    long start = clock();
    cout << "Doing pigeonhole process..." << endl;
    //filterReads(readsFilename, q, m, j, k);
    //parallelizeFilterReads(readsFilename, q, m, j, k);
    //selectingSeeds(readsFilename, q, m, j, k);
    selectingSeeds(readsFilename, windowLength, q);
    long end = clock();
    cout << "Done!" << endl;

    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by the pigeonhole process is : " << fixed
         << time_taken << " sec" << endl;

    return 0;
}


//    /**
//     * Direct Addressing
//     */
//    double dirTableSize = pow(4, q) + 1;
//    double posTableSize = refGenome.size() - q + 1;
//
//    for (int i = 0; i < dirTableSize; i++) {
//        dirTable.push_back(0);
//    }
//
//    for (int i = 0; i < posTableSize; i++) {
//        posTable.push_back(0);
//    }
//
//    buildTablesDirect(refGenome, q, dirTableSize, posTableSize, dirTable, posTable);
//
//
//    /**
//     * Open Addressing
//     */
//    double loadFactor = 0.8;
//    double codeTableSize = floor(( pow(loadFactor, -1)) * refGenome.size());
//    double dirTableSize = codeTableSize + 1;
//    double posTableSize = refGenome.size() - q + 1;
//    unsigned long int shiftedValue = ((unsigned long int)1 << (q * 2));
//
//    for (int i = 0; i < codeTableSize; i++) {
//        codeTable.push_back(-1);
//    }
//
//    for (int i = 0; i < dirTableSize; i++) {
//        dirTable.push_back(0);
//    }
//
//    for (int i = 0; i < posTableSize; i++) {
//        posTable.push_back(0);
//    }
//
//    buildTablesOpen(refGenome, q, shiftedValue, codeTableSize, dirTableSize, posTableSize, codeTable, dirTable, posTable);