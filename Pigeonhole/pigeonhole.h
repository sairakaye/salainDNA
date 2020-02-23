//
// Created by saimanalili on 30/12/2019.
//

#ifndef MP_PIGEONHOLE_PIGEONHOLE_H
#define MP_PIGEONHOLE_PIGEONHOLE_H

#include "common.h"

//typedef struct {
//    int startPos;
//    string read;
//    string readFromGenome;
//} AlignedReads;

extern int numberSeeds;
extern int numberAcceptedSeeds;
extern unsigned int numberOfFoundSeeds;
extern int numberReads;
extern int numberAcceptedReads;

string reverseComplement(string read);
void verification(vector<unsigned long long int>& forwardFound, vector<unsigned long long int>& reverseFound, vector<unsigned long long int>& exactFound);

#endif //MP_PIGEONHOLE_PIGEONHOLE_H