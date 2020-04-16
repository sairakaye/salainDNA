//
// Created by saimanalili on 25/02/2020.
//

#include "minimizers.h"

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

multimap<vector<unsigned long long int>, unsigned long long int>  generateMinimizers(string stringDNA, string mainName, unsigned int q, unsigned int w, unsigned int m)
{
    auto indexing_time_start = chrono::high_resolution_clock::now();

    ofstream outfile;
    string filename ("min_" + mainName + "_" + to_string(q) + ".txt");
    outfile.open(filename.c_str(), ios::out);

    unsigned long int noOfKmers = m - q + 1;
    unsigned long int windowSize = w + q - 1;

    cout << "Indexing start. \n";

    map<unsigned long long int, vector<unsigned long long int>> minimizers; // <hash key of minimizer, positions in ref>
    unsigned long long int mask = pow(4, q);

    for (int i = 0; i < (m - windowSize + 1); i++) {
        string finalMin;
        int finalMinIndex;

        string windowStr = stringDNA.substr(i, windowSize);
        unsigned long long int minm1 = 1L << (2 * q + 1);

        for (int j = 0; j < (windowSize - q + 1); j++){
            string sMinimizer = windowStr.substr(j, q);
            unsigned long long int hashValue = extractRanking(sMinimizer);
            unsigned long long int tempMinHash = inthash_64(hashValue, mask - 1);

            if (tempMinHash < minm1) {
                minm1 = tempMinHash;
                finalMin = sMinimizer;
                finalMinIndex = j + i + 1;
            }
        }

        unsigned long long int rankHashValue = extractRanking(finalMin);
        unsigned long long int finalMinHash = inthash_64(rankHashValue, mask - 1);
        minimizers[finalMinHash].push_back(i);
//      cout << "[" <<  windowStr << " " << i << "] ["  << finalMin << " " << (finalMinIndex - 1) << " " << finalMinHash << "]" << endl;
    }

    for (int i = 0; i < (w - 1); i++) {
        string finalMin;
        int finalMinIndex;

        string windowStr = stringDNA.substr(stringDNA.length() - (w + i), windowSize);
        unsigned long long int minm1 = 1L << (2 * q + 1);

        for (int j = 0; j < (windowStr.length() - q + 1); j++){
            string sMinimizer = windowStr.substr(j, q);
            unsigned long long int hashValue = extractRanking(sMinimizer);
            unsigned long long int tempMinHash = inthash_64(hashValue, mask - 1);

            if (tempMinHash < minm1) {
                minm1 = tempMinHash;
                finalMin = sMinimizer;
                finalMinIndex = stringDNA.length() - (w + i) + j;
            }
        }

        unsigned long long int rankHashValue = extractRanking(finalMin);
        unsigned long long int finalMinHash = inthash_64(rankHashValue,  mask - 1);
        minimizers[finalMinHash].push_back(stringDNA.length() - (w + i));
//		cout << "[" << windowStr << " " << stringDNA.length() - (w + i) << "] ["  << finalMin << " " << finalMinIndex << " " << finalMinHash << "]" << endl;
    }

    cout << "Indexing done. \n";

    auto indexing_time_end = chrono::high_resolution_clock::now();
    duration<double, std::milli> indexing_duration = duration_cast<milliseconds>(indexing_time_end - indexing_time_start);
    cout << "Total indexing time: " << (indexing_duration.count()/1000) << " s \n\nWriting start.\n";

    auto writing_time_start = chrono::high_resolution_clock::now();

    multimap<vector<unsigned long long int>, unsigned long long int> finalMinimizers;
    for(auto const &kv : minimizers)
        finalMinimizers.insert(make_pair(kv.second, kv.first));

    for (auto outerItr = finalMinimizers.begin(); outerItr != finalMinimizers.end(); outerItr++) {
        outfile << outerItr->second << ": ";

        for (auto innerItr = outerItr->first.begin(); innerItr != outerItr->first.end(); innerItr++) {
            outfile << *innerItr << " ";
        }
        outfile << endl;
    }

    cout << "Writing done. \n";

    auto writing_time_end = chrono::high_resolution_clock::now();
    duration<double, std::milli> writing_duration = duration_cast<milliseconds>(writing_time_end - writing_time_start);
    cout << "Total writing time: " << (writing_duration.count()/1000) << " s \n\n";

    cout << "No. of seeds: " << noOfKmers;
    cout << "\nNo. of minimizers: " << finalMinimizers.size();
    cout << "\nImprovement factor: " << ((long double)(m - q + 1) / minimizers.size())<< '\n';

    outfile.close();

    return finalMinimizers;
}

void buildMinimizersIndexing(string& genome, string& mainName) {
    unsigned long long int q2 = q;
    unsigned long long int w = q2;

    cout << "Minimizers for q = " << q << " \n\n";
    multimap<vector<unsigned long long int>, unsigned long long int> minimizers = generateMinimizers(genome, mainName, q, w, m);
}