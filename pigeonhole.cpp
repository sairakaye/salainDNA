//
// Created by saimanalili on 25/02/2020.
//

#include "pigeonhole.h"

int j;
int allowableE;

string reverseComplement(string read) {
    string reverseRead = read;
    reverse(reverseRead.begin(), reverseRead.end());

    for (int i = 0; i < reverseRead.length(); i++) {
        switch (reverseRead[i]) {
            case 'A':
                reverseRead[i] = 'T';
                break;
            case 'C':
                reverseRead[i] = 'G';
                break;
            case 'G':
                reverseRead[i] = 'C';
                break;
            case 'T':
                reverseRead[i] = 'A';
                break;
        }
    }

    return reverseRead;
}

void searchingReadProcess() {
    j = m / q;
    allowableE = floor(e / j);

    bool isCheckForApproximate = true;

    if (allowableE == 0) {
        isCheckForApproximate = false;
    }

    int i;
    #pragma omp parallel for
    for (i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i]);
        unsigned int tempAcceptedSeeds = 0;
        unsigned int tempReverseSeeds = 0;
        bool isFound = false;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);
            } else {
                seed = forwardRead.substr(startPosition, q);
            }

            if (seed.size() == w || seed.size() == q) {
                vector<unsigned long long int> forward;

                searchingPosition(seed, forwardRead, mode, q, k, true, isCheckForApproximate, forward);

                if (forward.size() > 0) {
                    tempAcceptedSeeds++;
                    isFound = true;
                }
            } else {
                break;
            }
        }

        if (!isFound) {
            string reverseRead = reverseComplement(forwardRead);

            int k;
            for (k = 0; k < j; k++) {
                tempReverseSeeds++;
                int startPosition = k * q;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);
                } else {
                    seed = reverseRead.substr(startPosition, q);
                }

                if (seed.size() == w || seed.size() == q) {
                    vector<unsigned long long int> reverse;

                    searchingPosition(seed, reverseRead, mode, q, k, false, isCheckForApproximate, reverse);

                    if (reverse.size() > 0) {
                        tempAcceptedSeeds++;
                    }
                } else {
                    break;
                }
            }
        }

        #pragma omp atomic
        counter->numAcceptedSeeds = counter->numAcceptedSeeds + tempAcceptedSeeds;
        #pragma omp atomic
        counter->numSeeds = counter->numSeeds + tempReverseSeeds;
    }
}

void searchingReadFoundExitProcess() {
    j = m / q;
    allowableE = floor(e / j);

    bool isCheckForApproximate = true;

    if (allowableE == 0) {
        isCheckForApproximate = false;
    }

    int i;
    #pragma omp parallel for
    for (i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i]);
        unsigned int tempAcceptedSeeds = 0;
        unsigned int tempReverseSeeds = 0;
        bool isFound = false;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);
            } else {
                seed = forwardRead.substr(startPosition, q);
            }

            if (seed.size() == w || seed.size() == q) {
                vector<unsigned long long int> forward;

                searchingPosition(seed, forwardRead, mode, q, k, true, isCheckForApproximate, forward);

                if (forward.size() > 0) {
                    tempAcceptedSeeds++;
                    isFound = true;
                    break;
                }
            } else {
                break;
            }
        }

        if (!isFound) {
            string reverseRead = reverseComplement(forwardRead);

            int k;
            for (k = 0; k < j; k++) {
                tempReverseSeeds++;
                int startPosition = k * q;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);
                } else {
                    seed = reverseRead.substr(startPosition, q);
                }

                if (seed.size() == w || seed.size() == q) {
                    vector<unsigned long long int> reverse;

                    searchingPosition(seed, reverseRead, mode, q, k, false, isCheckForApproximate, reverse);

                    if (reverse.size() > 0) {
                        tempAcceptedSeeds++;
                    }
                } else {
                    break;
                }
            }
        }

       #pragma omp atomic
        counter->numAcceptedSeeds = counter->numAcceptedSeeds + tempAcceptedSeeds;
       #pragma omp atomic
        counter->numSeeds = counter->numSeeds + tempReverseSeeds;
    }
}

