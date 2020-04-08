//
// Created by saimanalili on 25/02/2020.
//

#include "minimizer.h"

vector<pair<string, int> > alphabetRef = { {"A", 0}, {"C",1}, {"G",2}, {"T", 3} };

unsigned long long int extractRanking(string kMer) {
    string binary;
    int rankValue;

    for (int i = 0; i<kMer.length(); i++){
        for (int j = 0; j< alphabetRef.size(); j++){
            if (kMer.at(i) + string() == alphabetRef.at(j).first){
                rankValue = alphabetRef.at(j).second;
                binary.append(bitset<2>(rankValue).to_string());
            }
        }
    }

    unsigned long long int decimal = strtoull(binary.c_str(), nullptr, 2);
    return decimal;
}

uint64_t inthash_64(uint64_t key, uint64_t mask) {
    key = (~key + (key << 21)) & mask;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

unsigned long long int getMinimizerRank(string windowSeed, int q, int windowSize) {
    unsigned long long int mask = pow(4, q);
    string finalMin;

    unsigned long long int minm1 = 1L << (2 * q + 1);

    for (unsigned int j = 0; j < (windowSize - q + 1); j++){
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