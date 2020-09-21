//
// Created by saimanalili on 25/02/2020.
//

#include "minimizers.h"

unsigned long long int getMinimizerRank(string windowSeed, int q, int windowSize) {
    unsigned long long int mask = pow(4, q);
    string finalMin;

    unsigned long long int minm1 = 1ULL << (2 * q + 1);

    for (unsigned int j = 0; j < (windowSize - q + 1); j++) {
        string sMinimizer = windowSeed.substr(j, q);
        unsigned long long int hashValue = extractRanking(sMinimizer);
        unsigned long long int tempMinHash = inthash_64(hashValue, mask - 1);

        if (tempMinHash < minm1) {
            minm1 = tempMinHash;
            finalMin = sMinimizer;
        }
    }

    unsigned long long int rankHashValue = extractRanking(finalMin);
    unsigned long long int finalMinHash = inthash_64(rankHashValue, mask - 1);

    return finalMinHash;
}

unsigned long long int getMinimizerRankWithoutWindow(string windowSeed, int q) {
    unsigned long long int mask = pow(4, q);

    unsigned long long int rankHashValue = extractRanking(windowSeed);
    unsigned long long int finalMinHash = inthash_64(rankHashValue, mask - 1);

    return finalMinHash;
}

multimap<vector<unsigned long long int>, unsigned long long int> generateMinimizers(string genome, string genomeFileName, unsigned int q, unsigned int w, unsigned int n)
{
    ofstream outfile;
    string filename ("min_" + genomeFileName + "_" + to_string(q) + ".txt");
    outfile.open(filename.c_str(), ios::out);

    unsigned long int windowSize = w + q - 1;

    map<unsigned long long int, vector<unsigned long long int>> initMinimizers;
    unsigned long long int mask = pow(4, q);

    unsigned int i;

    #pragma omp parallel for
    for (i = 0; i < n - windowSize + 1; i++) {
        string minimizer;

        string window = genome.substr(i, windowSize);
        unsigned long long int min = 1ULL << (2 * q + 1);

        for (unsigned int j = 0; j < windowSize - q + 1; j++){
            string qGram = window.substr(j, q);
            unsigned long long int rank = extractRanking(qGram);
            unsigned long long int hashValue = inthash_64(rank, mask - 1);

            if (hashValue < min) {
                min = hashValue;
                minimizer = qGram;
            }
        }

        #pragma omp critical
        {
            unsigned long long int finalRank = extractRanking(minimizer);
            unsigned long long int finalHashValue = inthash_64(finalRank, mask - 1);
            initMinimizers[finalHashValue].push_back(i);
        }
    }

    #pragma omp parallel for
    for (i = 0; i < w - 1; i++) {
        string minimizer;

        string window = genome.substr(n - (w + i), windowSize);
        unsigned long long int min = 1ULL << (2 * q + 1);

        for (unsigned int j = 0; j < window.length() - q + 1; j++){
            string qGram = window.substr(j, q);
            unsigned long long int rank = extractRanking(qGram);
            unsigned long long int hashValue = inthash_64(rank, mask - 1);

            if (hashValue < min) {
                min = hashValue;
                minimizer = qGram;
            }
        }

        #pragma omp critical
        {
            unsigned long long int finalRank = extractRanking(minimizer);
            unsigned long long int finalHashValue = inthash_64(finalRank, mask - 1);
            initMinimizers[finalHashValue].push_back(n - (w + i));
        }
    }

    multimap<vector<unsigned long long int>, unsigned long long int> finalMinimizers;
    for(auto const &kv : initMinimizers)
        finalMinimizers.insert(make_pair(kv.second, kv.first));

    for (auto outerItr = finalMinimizers.begin(); outerItr != finalMinimizers.end(); outerItr++) {
        outfile << outerItr->second << ": ";

        for (auto innerItr = outerItr->first.begin(); innerItr != outerItr->first.end(); innerItr++) {
            outfile << *innerItr << " ";
            minimizers[outerItr->second].push_back(*innerItr);
        }
        outfile << endl;
    }

    outfile.close();

    return finalMinimizers;
}

void buildMinimizersIndexingFile(string& genome, string& genomeFileName) {
    unsigned long long int qIndexing = (unsigned long long int)q;
    unsigned long long int wIndexing = qIndexing;

    cout << "Minimizer-based for q = " << qIndexing << endl << endl;
    multimap<vector<unsigned long long int>, unsigned long long int> minimizers = generateMinimizers(genome, genomeFileName, qIndexing, wIndexing, genome.length());
}