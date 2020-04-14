//
// Created by saimanalili on 13/04/2020.
//

#include "indexing.h"

int occupiedSpaces = 0;
int collisions = 0;

vector<string> qGrams;

void generateQGrams(string prefix, int k)
{
    char charactersDNA[] = {'A', 'C', 'G', 'T'};

    if (k == 0) {
        qGrams.push_back(prefix);
        return;
    }

    for (int i = 0; i < 4; i++) {
        string newPrefix;
        newPrefix = prefix + charactersDNA[i];
        generateQGrams(newPrefix, k - 1);
    }
}

void buildDirectAddressingTables(string stringDNA, string mainName, int m, int q, unsigned long long int dirTableSize, unsigned long long int posTableSize)
{
    auto indexing_time_start = chrono::high_resolution_clock::now();

    ofstream outfile;
    outfile.open("dir_" + mainName + "_" + to_string(q) + ".txt", ios::out);

    vector<unsigned long> dirTable(dirTableSize);
    vector<unsigned long> posTable(posTableSize);

    fill(dirTable.begin(), dirTable.end(), 0);
    fill(posTable.begin(), posTable.end(), 0);

    cout << "Indexing start. \n";

    generateQGrams("", q);

//	vector<string> qGramsFinal(qGrams.size()); //
//	for (int i = 0; i < qGrams.size(); i++) { //
//		int kMerIndexInGenome = extractRanking(qGrams.at(i)); //
//
//		int hashValue = inthash_64(kMerIndexInGenome, dirTableSize - 2); //
//		qGramsFinal.at(hashValue) = qGrams.at(i); //
//	}

    for (int i = 0; i < (m - (q - 1)); i++) {
        unsigned int tempIndex = extractRanking(stringDNA.substr(i, q));

        cout << stringDNA.substr(i, q)  << ' ' << tempIndex << endl;
        dirTable.at(tempIndex)++;
    }

    for (int i = 1; i < dirTable.size(); i++) {
        dirTable.at(i) += dirTable.at(i-1);
    }

    int index = 0;
    for (int i = m - q; i >= 0; i--) {
        unsigned long int tempIndex = extractRanking(stringDNA.substr(i, q));
        cout << stringDNA.substr(i, q)  << ' ' << tempIndex << endl;
        dirTable.at(tempIndex)--;

        index = dirTable.at(tempIndex);
        posTable.at(index) = i;
    }

    cout << "Indexing done. \n";

    auto indexing_time_end = chrono::high_resolution_clock::now();
    duration<double, std::milli> indexing_duration = duration_cast<milliseconds>(indexing_time_end - indexing_time_start);
    cout << "Total indexing time: " << (indexing_duration.count()/1000) << " s \n\nWriting start.\n";

    auto writing_time_start = chrono::high_resolution_clock::now();

    outfile << "dir" << "\n";
//	cout << "dir Table [of size " << dirTableSize << "] : ";
    for (int i = 0; i < dirTable.size(); i++) {
//		cout << dirTable.at(i) << ' ';
        outfile << dirTable.at(i) << "\n";
    }
    outfile << "\n" << "pos" << "\n";
//	cout << "\npos Table [of size " << posTableSize << "] : ";
    for (int i = 0; i < posTable.size(); i++) {
//		cout << posTable.at(i) << ' ';
        outfile << posTable.at(i) << "\n";
    }
//    outfile << "\n";

    cout << "Writing done. \n";

    auto writing_time_end = chrono::high_resolution_clock::now();
    duration<double, std::milli> writing_duration = duration_cast<milliseconds>(writing_time_end - writing_time_start);
    cout << "Total writing time: " << (writing_duration.count()/1000) << " s \n\n";

    outfile.close();
}

void buildOpenAddressingTables(string stringDNA, string mainName, int m, int q, unsigned int codeTableSize, unsigned int dirTableSize, int posTableSize) {
    auto indexing_time_start = chrono::high_resolution_clock::now();

    ofstream outfile;
    outfile.open("open_" + mainName + "_" + to_string(q) + ".txt", ios::out);

    vector<long long> codeTable(codeTableSize);
    vector<unsigned long> dirTable(dirTableSize);
    vector<unsigned long> posTable(posTableSize);

    fill(codeTable.begin(), codeTable.end(), -1);
    fill(dirTable.begin(), dirTable.end(), 0);
    fill(posTable.begin(), posTable.end(), 0);

    vector<int> codeTableIndices; // store the indices here

    cout << "Indexing start. \n";

    for (int i = 0; i < (m-q+1); i++) {
        unsigned long int kMerIndexInGenome = extractRanking(stringDNA.substr(i, q));
        unsigned int hashValue = kMerIndexInGenome % (int)codeTableSize;
//		int hashValue = inthash_64(kMerIndexInGenome, shiftedValue - 1);

        if (codeTable.at(hashValue) != -1) {
            collisions++;
            unsigned int j = hashValue;
            unsigned int k = 1;

            while (codeTable.at(j) != -1 && codeTable.at(j) != kMerIndexInGenome) {
                j = (j + (k * k)) % (int)codeTableSize;
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
        unsigned long int kMerIndexInGenome = extractRanking(stringDNA.substr(i, q));
        dirTable.at(codeTableIndices.at(i))++;
    }

    for (int i = 1; i < dirTable.size(); i++) {
        dirTable.at(i) += dirTable.at(i-1);
    }

    unsigned int index = 0;
    for (int i = stringDNA.length() - q; i >= 0; i--) {
        unsigned long int kMerIndexInGenome = extractRanking(stringDNA.substr(i, q));
        dirTable.at(codeTableIndices.at(i))--;
        index = dirTable.at(codeTableIndices.at(i));
        posTable.at(index) = i;
    }

    cout << "Indexing done. \n";

    auto indexing_time_end = chrono::high_resolution_clock::now();
    duration<double, std::milli> indexing_duration = duration_cast<milliseconds>(indexing_time_end - indexing_time_start);
    cout << "Total indexing time: " << (indexing_duration.count()/1000) << " s \n\nWriting start.\n";

    auto writing_time_start = chrono::high_resolution_clock::now();

    outfile << "code" << endl;
//	cout << "code Table [of size " << codeTableSize << "] : ";
    for (int i = 0; i < codeTable.size(); i++) {
//		cout << codeTable.at(i) << ' ';

        if (codeTable.at(i) != -1) {
            outfile << (unsigned long int)codeTable.at(i) << " " << i << endl;
        } else {
            outfile << codeTable.at(i) << " " << i << endl;
        }

//        if (codeTable.at(i) != (-1)) {
//            outfile << i << " " << (unsigned int)codeTable.at(i) << endl;
//        } else {
//            outfile << i << " " << codeTable.at(i) << endl;
//        }
    }
    outfile << endl << "dir" << endl;
//	cout << "\ndir Table [of size " << dirTableSize << "] : ";
    for (int i = 0; i < dirTable.size(); i++) {
//		cout << dirTable.at(i) << ' ';
        outfile << dirTable.at(i) << endl;
    }
    outfile << endl << "pos" << endl;
//	cout << "\npos Table [of size " << posTableSize << "] : ";
    for (int i = 0; i < posTable.size(); i++) {
//		cout << posTable.at(i) << ' ';
        outfile << posTable.at(i) << endl;
    }
//	cout << endl;
    outfile << endl;

    cout << "Writing done. \n";

    auto writing_time_end = chrono::high_resolution_clock::now();
    duration<double, std::milli> writing_duration = duration_cast<milliseconds>(writing_time_end - writing_time_start);
    cout << "Total writing time: " << (writing_duration.count()/1000) << " s \n\n";

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


void buildDirectAddressingIndexing(string& genome, string& mainName) {
    unsigned long long int dirTableSize = pow(4, q) + 1;
    int posTableSize = genome.size() - q + 1;

    cout << "Direct Addressing for q = " << q << " \n\n";
    buildDirectAddressingTables(genome, mainName, genome.length(), q, dirTableSize, posTableSize);

//		auto total_time_end = chrono::high_resolution_clock::now();
//		duration<double, std::milli> total_duration = duration_cast<milliseconds>(total_time_end - total_time_start);
//		cout << "\nTotal execution time: " << total_duration.count() << " ms \n";

    cout << "Done.\n";
}

void buildOpenAddressingIndexing(string& genome, string& mainName) {
    double loadFactor = 0.8;

    unsigned int codeTableSize = floor(( pow(loadFactor, -1)) * genome.size());
    unsigned int dirTableSize = codeTableSize + 1;
    unsigned int posTableSize = genome.size() - q + 1;

//		cout << codeTableSize << ' '<< dirTableSize << ' ' << ' ' << posTableSize;

    cout << "Open Addressing for q = " << q << "\n\n";
    buildOpenAddressingTables(genome, mainName, genome.length(), q, codeTableSize, dirTableSize, posTableSize);

//		auto total_time_end = chrono::high_resolution_clock::now();
//		duration<double, std::milli> total_duration = duration_cast<milliseconds>(total_time_end - total_time_start);
//		cout << "\nTotal execution time: " << (total_duration.count()/1000) << " s \n";

//		double percentageOccupied = ((double)occupiedSpaces / codeTableSize) * 100;
//		double percentageCollided = ((double)collisions / occupiedSpaces) * 100;
//		cout << "\n" << occupiedSpaces << " of " << codeTableSize << " places taken (" << percentageOccupied << "%)";
//		cout << "\n" << collisions << " of " << occupiedSpaces << " collisions ("  << percentageCollided << "%)" << endl;

    cout << "Done.\n";
}

void buildMinimizersIndexing(string& genome, string& mainName) {
    unsigned long long int q2 = q;
    unsigned long long int w = q2;
    unsigned long long int m = genome.length();

    cout << "Minimizers for q = " << q << " \n\n";
    multimap<vector<unsigned long long int>, unsigned long long int> minimizers = generateMinimizers(genome, mainName, q, w, m);
}