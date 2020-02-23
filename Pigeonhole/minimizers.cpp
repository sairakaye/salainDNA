//
// Created by saimanalili on 08/01/2020.
//

#include "minimizers.h"

vector<pair<string, int> > alphabetRef = { {"A", 0}, {"C",1}, {"G",2}, {"T", 3} };

unsigned long long int extractRanking(string kMer)
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

    unsigned long long int decimal = std::strtoul (binary.c_str(), nullptr, 2);
    return decimal;
}

uint64_t inthash_64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

unsigned int getMinimizerRank(string windowSeed, int q, int windowSize) {
    double mask = double(pow(4, q));
    string finalMin;

    unsigned long long int minm1 = 1L << (2 * q + 1);

    for (unsigned int j = 0; j < (windowSize - q + 1); j++){
        string sMinimizer = windowSeed.substr(j, q);
        unsigned long long int hashValue = extractRanking(sMinimizer);
        int tempMinHash = inthash_64(hashValue, mask - 1);

        if (tempMinHash < minm1) {
            minm1 = tempMinHash;
            finalMin = sMinimizer;
        }
    }

    unsigned long long int rankHashValue = extractRanking(finalMin);
    unsigned int finalMinHash = inthash_64(rankHashValue, mask - 1);

    return finalMinHash;
}

/*
uint64_t inthash_64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}
multimap<vector<int>, int> generateMinimizers(string stringDNA, int q, int w, int m) {
//	ofstream outfile;
//	outfile.open("min_1m.txt", ios::out);
    unsigned long long int noOfKmers = m - q + 1;
    unsigned long long int windowSize = w + q - 1;
    cout << "\n\nFormat: \n[seed, seed index]  [minimizer, minimizer index, int hash key] \n\n";
    map<int, vector<int>> minimizers; // <hash key of minimizer, positions in ref>
    double mask = pow(4, q);
    cout << "> Minimizers: " << endl;
    for (int i = 0; i < (m - windowSize + 1); i++) {
        string finalMin;
        int finalMinIndex;
        string windowStr = stringDNA.substr(i, windowSize);
        unsigned long long int minm1 = 1 << (2 * q + 1);
        for (int j = 0; j < (windowSize - q + 1); j++){
            string sMinimizer = windowStr.substr(j, q);
            int hashValue = extractRanking(sMinimizer);
            int tempMinHash = inthash_64(hashValue, mask - 1);
            if (tempMinHash < minm1) {
                minm1 = tempMinHash;
                finalMin = sMinimizer;
                finalMinIndex = j + i + 1;
            }
        }
        int rankHashValue = extractRanking(finalMin);
        int finalMinHash = inthash_64(rankHashValue, mask - 1);
        minimizers[finalMinHash].push_back(i);
        cout << "[" <<  windowStr << " " << i << "] ["  << finalMin << " " << (finalMinIndex - 1) << " " << finalMinHash << "]" << endl;
    }
    cout << "> Rear-end Minimizers: " << endl;
    for (int i = 0; i < (w - 1); i++) {
        string finalMin;
        int finalMinIndex;
        string windowStr = stringDNA.substr(stringDNA.length() - (w + i));
        unsigned long long int minm1 = 1 << (2 * q + 1);
        for (int j = 0; j < (windowStr.length() - q + 1); j++){
            string sMinimizer = windowStr.substr(j, q);
            int hashValue = extractRanking(sMinimizer);
            int tempMinHash = inthash_64(hashValue, mask - 1);
            if (tempMinHash < minm1) {
                minm1 = tempMinHash;
                finalMin = sMinimizer;
                finalMinIndex = stringDNA.length() - (w + i) + j;
            }
        }
        int rankHashValue = extractRanking(finalMin);
        int finalMinHash = inthash_64(rankHashValue,  mask - 1);
        minimizers[finalMinHash].push_back(stringDNA.length() - (w + i));
        cout << "[" << windowStr << " " << stringDNA.length() - (w + i) << "] ["  << finalMin << " " << finalMinIndex << " " << finalMinHash << "]" << endl;
    }
    multimap<vector<int>, int> finalMinimizers;
    for(auto const &kv : minimizers)
        finalMinimizers.insert(make_pair(kv.second, kv.first));
    cout << "\nFinal minimizers: (int hash key: seed indices)" << endl;
    for (auto outerItr = finalMinimizers.begin(); outerItr != finalMinimizers.end(); outerItr++) {
        cout << outerItr->second << ": ";
//	  outfile << outerItr->second << ": ";
        for (auto innerItr = outerItr->first.begin(); innerItr != outerItr->first.end(); innerItr++) {
            cout << *innerItr << " ";
//		  outfile << *innerItr << " ";
        }
        cout << endl;
//	  outfile << endl;
    }
    cout << "\nNo. of seeds: " << noOfKmers;
    cout << "\nNo. of minimizers: " << finalMinimizers.size();
    cout << "\nImprovement factor: " << ((long double)(m - q + 1) / minimizers.size());
//	outfile.close();
    return finalMinimizers;
}*/