//
// Created by saimanalili on 25/02/2020.
//

#include "pigeonhole.h"

int j;
int allowableE;

void searchingReadProcess() {
    j = m / q;
    allowableE = floor(e / j);

    bool isExactMatching = false;

    if (allowableE == 0) {
        isExactMatching = true;
    }

    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i]);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);

                if (seed.size() < w) {
                    startPosition += seed.size() - w;
                    seed = forwardRead.substr(startPosition, w);
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.size() < q) {
                    startPosition += seed.size() - q;
                    seed = forwardRead.substr(startPosition, q);
                }
            }

            vector<unsigned long long int> forward;

            searchingPosition(seed, forwardRead, mode, q, k, true, isExactMatching, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());
            }
        }

        if (isExactMatching && tempAcceptedSeeds >= (j - e)) {
            #pragma omp critical
            {
                possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
            };
        } else if (isExactMatching && tempAcceptedSeeds == 0 && isReverseAccepted) {
            tempAcceptedSeeds = 0;
            string reverseRead = reverseComplement(forwardRead);

            int k;
            for (k = 0; k < j; k++) {
                int startPosition = k * q;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);

                    if (seed.size() < w) {
                        startPosition += seed.size() - w;
                        seed = forwardRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.size() < q) {
                        startPosition += seed.size() - q;
                        seed = forwardRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                searchingPosition(seed, reverseRead, mode, q, k, false, isExactMatching, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());
                }
            }

            if (tempAcceptedSeeds >= (j - e)) {
                #pragma omp critical
                {
                    possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
                };
            }
        } else if (!isExactMatching && tempAcceptedSeeds == j) {
            #pragma omp critical
            {
                possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
            };

            #pragma omp atomic
            numAcceptedReads++;
        } else if (!isExactMatching && tempAcceptedSeeds == 0 && isReverseAccepted) {
            tempAcceptedSeeds = 0;
            string reverseRead = reverseComplement(forwardRead);

            int k;
            for (k = 0; k < j; k++) {
                int startPosition = k * q;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);

                    if (seed.size() < w) {
                        startPosition += seed.size() - w;
                        seed = forwardRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.size() < q) {
                        startPosition += seed.size() - q;
                        seed = forwardRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                searchingPosition(seed, reverseRead, mode, q, k, false, isExactMatching, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());
                }
            }

            if (tempAcceptedSeeds == j) {
                #pragma omp critical
                {
                    possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
                };
            }
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

/*
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

                if (seed.size() < w) {
                    startPosition = startPosition + seed.size() - w;
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.size() < q) {
                    startPosition = startPosition + seed.size() - q;
                }
            }

            vector<unsigned long long int> forward;

            searchingPosition(seed, forwardRead, mode, q, k, true, isCheckForApproximate, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                isFound = true;
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
        numAcceptedSeeds += tempAcceptedSeeds;
        #pragma omp atomic
        numSeeds += tempReverseSeeds;
    }
}
*/
