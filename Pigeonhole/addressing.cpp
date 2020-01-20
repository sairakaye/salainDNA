//============================================================================
// Name        : DNADirectAddressing.cpp
// Author      : Candace Claire Mercado
// Version     :
// Copyright   : Your copyright notice
// Description : DNA Direct Addressing Q-gram Index in C++, Ansi-style
//============================================================================

#include "addressing.h"

vector<string> qGrams;

bool isFound (vector<string> vec, string ele){
    if (find(vec.begin(), vec.end(), ele) != vec.end())
        return true;
    else return false;
}

int indexOf (vector<string> vec, string ele)
{
    int indexFoundAt = (find(vec.begin(), vec.end(), ele) - vec.begin());

    if(indexFoundAt >= vec.size()) {
        return -1;
    } else {
        return indexFoundAt;
    }
}

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

void buildTablesDirect(string stringDNA, int q, double dirTableSize, double posTableSize, vector<int> dirTable, vector<int> posTable)
{
//	vector<string> kMersInGenome;
//	for (int i = 0; i < stringDNA.length() - (q - 1); i++) {
//		if (!isFound(kMersInGenome, stringDNA.substr(i, q)))
//			kMersInGenome.push_back(stringDNA.substr(i, q));
//	}

    generateQGrams("", q);

//	vector<string> qGramsFinal(qGrams.size()); //
//	for (int i = 0; i < qGrams.size(); i++) { //
//		int kMerIndexInGenome = extractRanking(qGrams.at(i)); //
//
//		int hashValue = inthash_64(kMerIndexInGenome, dirTableSize - 2); //
//		qGramsFinal.at(hashValue) = qGrams.at(i); //
//	}

    for (int i = 0; i < stringDNA.length() - (q - 1); i++) {
        int tempIndex = indexOf(qGrams, stringDNA.substr(i, q));
        dirTable.at(tempIndex)++;
    }

    for (int i = 1; i < dirTable.size(); i++) {
        dirTable.at(i) += dirTable.at(i-1);
    }

    int index = 0;
    for (int i = stringDNA.length() - q; i >= 0; i--) {
        int tempIndex = indexOf(qGrams, stringDNA.substr(i, q));
        dirTable.at(tempIndex)--;

        index = dirTable.at(tempIndex);
        posTable.at(index) = i;
    }

    cout << "dir Table [of size " << dirTableSize << "] : ";
    for (auto i : dirTable)
        cout << i << ' ';
    cout << "\npos Table [of size " << posTableSize << "] : ";
    for (auto i : posTable)
        cout << i << ' ';
}

/*int main(int argc, char **argv) {
    vector<string> genomeStrings;
    genomeStrings = initGenomeFile(genomeStrings);

    unsigned long int q = 12;

    if (genomeStrings[0].at(0) == '>') {
        string stringDNA;
        for (int i = 1; i < genomeStrings.size(); i++)
            stringDNA.append(genomeStrings[i]);


        unsigned long int shiftedValue = ((unsigned long int)1 << (q * 2));
//		unsigned long int shiftedValue = pow((unsigned long int)4, q);
        cout << "Number of possible q-grams to be generated: " << shiftedValue << "\n";

        double dirTableSize = pow(4, q) + 1;
        double posTableSize = stringDNA.size() - q + 1;

        cout << "Direct Addressing for q = " << q << "\n\n";
        buildTablesDirect(stringDNA, q, dirTableSize, posTableSize);

//		double percentageOccupied = ((double)occupiedSpaces / codeTableSize) * 100;
//		double percentageCollided = ((double)collisions / occupiedSpaces) * 100;
//
//		cout << "\n\n" << occupiedSpaces << " of " << codeTableSize << " places taken (" << percentageOccupied << "%)";
//		cout << "\n" << collisions << " of " << occupiedSpaces << " collisions ("  << percentageCollided << "%)";
    }
}*/

int indexOf (vector<int> vec, int ele){
    int indexFoundAt = (find(vec.begin(), vec.end(), ele) - vec.begin());
    if(indexFoundAt >= vec.size()) {
        return -1;
    } else {
        return indexFoundAt;
    }
}

int count (vector<string> vec, string ele){
    return count(vec.begin(), vec.end(), ele);
}

int count2 (vector<int> vec, int ele){
    return count(vec.begin(), vec.end(), ele);
}

void buildTablesOpen(string stringDNA, int q, int shiftedValue, double codeTableSize, double dirTableSize, double posTableSize, vector<int> codeTable, vector<int> dirTable, vector<int> posTable) {
    vector<string> kMersInGenome;
    for (int i = 0; i < stringDNA.length() - q + 1; i++) {
        if (!isFound(kMersInGenome, stringDNA.substr(i, q)))
            kMersInGenome.push_back(stringDNA.substr(i, q));
    }

    for (int i = 0; i < kMersInGenome.size(); i++) {
        int kMerIndexInGenome = extractRanking(kMersInGenome.at(i));

        int hashValue = inthash_64(kMerIndexInGenome, shiftedValue - 1);
//		cout << hashValue << " is the hash value for q-gram " << kMersInGenome.at(i) << "\n";

        if (codeTable.at(hashValue) == -1) {
            codeTable.at(hashValue) = kMerIndexInGenome;
        }
    }

    for (int i = 0; i < stringDNA.length() - q + 1; i++) {
        if (count(kMersInGenome, stringDNA.substr(i, q)) > 0){
            int kMerIndexInGenome = extractRanking(stringDNA.substr(i, q));
            dirTable.at(indexOf(codeTable, kMerIndexInGenome))++;
        }
    }

    for (int i = 1; i < dirTable.size(); i++) {
        dirTable.at(i) += dirTable.at(i-1);
    }

    int index = 0;
    for (int i = stringDNA.length() - q; i >= 0; i--) {
        int kMerIndexInGenome = extractRanking(stringDNA.substr(i, q));
        if (count2(codeTable, kMerIndexInGenome ) > 0) {
            dirTable.at(indexOf(codeTable, kMerIndexInGenome))--;
            index = dirTable.at(indexOf(codeTable, kMerIndexInGenome));
            posTable.at(index) = i;
        }
    }

    cout << "code Table [of size " << codeTableSize << "] : ";
    for (auto i : codeTable)
        cout << i << ' ';
    cout << "\ndir Table [of size " << dirTableSize << "] : ";
    for (auto i : dirTable)
        cout << i << ' ';
    cout << "\npos Table [of size " << posTableSize << "] : ";
    for (auto i : posTable)
        cout << i << ' ';
}


/*
int main(int argc, char **argv) {
    vector<string> genomeStrings;
    genomeStrings = initGenomeFile(genomeStrings);

    unsigned long int q = 8;
    double loadFactor = 0.8;

    if (genomeStrings[0].at(0) == '>') {
        string stringDNA;
        for (int i = 1; i < genomeStrings.size(); i++)
            stringDNA.append(genomeStrings[i]);

        unsigned long int shiftedValue = ((unsigned long int)1 << (q * 2));
//		unsigned long int shiftedValue = pow((unsigned long int)4, q);
        cout << "Number of possible q-grams to be generated: " << shiftedValue << "\n";

        double codeTableSize = floor(( pow(loadFactor, -1)) * stringDNA.size());
//		double codeTableSize =  stringDNA.size()+1;
        double dirTableSize = codeTableSize + 1;
        double posTableSize = stringDNA.size() - q + 1;

        cout << "Open Addressing for q = " << q << "\n\n";
        buildTablesOpen(stringDNA, q, shiftedValue, codeTableSize, dirTableSize, posTableSize);

        double percentageOccupied = ((double)occupiedSpaces / codeTableSize) * 100;
        double percentageCollided = ((double)collisions / occupiedSpaces) * 100;

        cout << "\n\n" << occupiedSpaces << " of " << codeTableSize << " places taken (" << percentageOccupied << "%)";
        cout << "\n" << collisions << " of " << occupiedSpaces << " collisions ("  << percentageCollided << "%)";
    }
}*/
