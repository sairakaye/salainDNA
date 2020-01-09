//
// Created by saimanalili on 29/12/2019.
//
#include <cstring>
#include "indexing.h"
#include "file_manager.h"

using namespace std;

string getCodedKmer(string kmer) {
    for (long i = 0; i < kmer.length(); i++) {
        switch (kmer[i]) {
            case 'A':
                kmer[i] = '0';
                break;

            case 'C':
                kmer[i] = '1';
                break;

            case 'G':
                kmer[i] = '2';
                break;

            case 'T':
                kmer[i] = '3';
                break;
        }
    }

    return kmer;
}

long getRanking(string kmer, int k) {
    long rank = 0;

    string translatedKmer = getCodedKmer(kmer);

    for (int i = 0; i < translatedKmer.length(); i++) {
        rank += (long)(translatedKmer[i] - '0') * pow(4, k - 1 - i);
    }

    return rank;
}

int findIndex(Code codeTable[], int codeTableLength, long element) {
    int index = -1;

    for (int i = 0; i < codeTableLength; i++) {
        if (codeTable[i].rank == element) {
            index = i;
            break;
        }
    }

    return index;
}

int quadraticProbing(long rank, HashIndexing *hashIndexing, int length) {
    int hashValue = rank % length;
    int codeTableIndex;

    if (hashIndexing->codeTable[hashValue].rank != -1 && findIndex(hashIndexing->codeTable, hashIndexing->codeTableLength, rank) < 0) {
        int j = hashValue;
        int k = 1;

        while (hashIndexing->codeTable[j].rank != -1) {
            j = (j + (k * k)) % length;
            k++;
        }

        hashIndexing->codeTable[j].rank = rank;
        codeTableIndex  = j;
    } else {
        hashIndexing->codeTable[hashValue].rank = rank;
        codeTableIndex = hashValue;
    }

    return codeTableIndex;
}

char *substring(char *s, int index, int n) {
    char *res = new char[n + 1];
    memcpy(res, s + index, n);
    res[n] = 0;
    return res;
}

void processIndexing(HashIndexing *hashIndexing, RefGenome *refGenome) {
    for (int i = 0; i < refGenome->length - (hashIndexing->k - 1); i++) {
        long rank = getRanking(substring(refGenome->genome, i, hashIndexing->k), hashIndexing->k);
        int codeTableIndex = quadraticProbing(rank, hashIndexing, hashIndexing->codeTableLength);
        //hashIndexing->codeTable[codeTableIndex].kmer = substring(refGenome->genome, i, hashIndexing->k);
        hashIndexing->dirTable[codeTableIndex]++;
    }

    for (int i = 1; i < hashIndexing->dirTableLength; i++) {
        hashIndexing->dirTable[i] += hashIndexing->dirTable[i-1];
    }

    for (int i = refGenome->length - hashIndexing->k; i >= 0; i--) {
        long rank = getRanking(substring(refGenome->genome, i, hashIndexing->k), hashIndexing->k);
        int index = findIndex(hashIndexing->codeTable, hashIndexing->codeTableLength, rank);

        if (index >= 0) {
            hashIndexing->dirTable[index]--;
            hashIndexing->posTable[hashIndexing->dirTable[index]] = i;
        }
    }
}