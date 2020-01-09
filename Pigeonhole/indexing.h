//
// Created by saimanalili on 29/12/2019.
//

#ifndef MP_PIGEONHOLE_INDEXING_H
#define MP_PIGEONHOLE_INDEXING_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <cstring>

#include "file_manager.h"

typedef struct {
    char *kmer = NULL;
    long rank;
} Code;


class HashIndexing {
public:
    Code *codeTable;
    int *dirTable;
    int *posTable;
    int k;
    int codeTableLength;
    int dirTableLength;
    int posTableLength;

    HashIndexing(int genomeLength, int k, float loadFactor) {
        this->k = k;
        codeTableLength = floor(pow(loadFactor, -1) * genomeLength);
        dirTableLength = codeTableLength + 1;
        posTableLength = genomeLength - k + 1;

        codeTable = (Code *)malloc(codeTableLength * sizeof(Code));

        if (codeTable == NULL) {
            std::cout << "Memory allocation failed." << std::endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < codeTableLength; i++) {
            (codeTable + i)->rank = -1;
        }

        dirTable = (int *)malloc(dirTableLength * sizeof(int));

        if (dirTable == NULL) {
            std::cout << "Memory allocation failed." << std::endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < dirTableLength; i++) {
            dirTable[i] = 0;
        }

        posTable = (int *)malloc(posTableLength * sizeof(int));

        if (posTable == NULL) {
            std::cout << "Memory allocation failed." << std::endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < posTableLength; i++) {
            posTable[i] = 0;
        }
    }

    ~HashIndexing() {
        if (this->codeTable != NULL) {
            free(this->codeTable);
        }

        if (this->dirTable != NULL) {
            free(this->dirTable);
        }

        if (this->posTable != NULL) {
            free(this->posTable);
        }
    }

    int findIndexOfKmer(const char *element) {
        // Optimize this. This is still O(n).

        int index = -1;

        for (int i = 0; i < this->codeTableLength; i++) {
            if (strcmp(this->codeTable[i].kmer, element) == 0) {
                index = i;
                break;
            }
        }

        return index;
    }
};

std::string getCodedKmer(std::string kmer);
long getRanking(std::string kmer, int k);
int findIndex(Code codeTable[], int codeTableLength, long element);
int quadraticProbing(long rank, HashIndexing *hashIndexing, int length);
char *substring(char *s, int index, int n);
void processIndexing(HashIndexing *hashIndexing, RefGenome *refGenome);

#endif //MP_PIGEONHOLE_INDEXING_H
