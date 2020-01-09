//
// Created by saimanalili on 08/01/2020.
//

//============================================================================
// Name        : DNAMinimizersIndex.cpp
// Author      : Candace Claire Mercado
// Version     :
// Copyright   : Your copyright notice
// Description : DNA Minimizer Implementation in C++, Ansi-style
//============================================================================

#include "minimizers.h"
using namespace std;

vector<pair<string, int> > alphabetRef = { {"A", 0}, {"C",1}, {"G",2}, {"T", 3} };

int extractRanking(string kMer)
{
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

    int decimal = stoi (binary, nullptr, 2);
    return decimal;
}

int generateMinimizerHash(string windowSeed, int q, int w) {
    unsigned long int windowSize = w + q - 1;

    string finalMin;

    unsigned long long int minm1 = 1 << (2 * q + 1);

    for (int j = 0; j < (windowSize - q + 1); j++) {
        string sMinimizer = windowSeed.substr(j, q);
        int hashValue = extractRanking(sMinimizer);

        if (hashValue < minm1) {
            minm1 = hashValue;
            finalMin = sMinimizer;
        }
    }

    int rankHashValue = extractRanking(finalMin);
    return rankHashValue;
}