//
// Created by saimanalili on 25/02/2020.
//

#include "pigeonhole.h"

int j;
int allowableE;

void exactSearchingForAll() {
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

            exactSearchingPosition(seed, forwardRead, mode, q, k, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());
            }
        }

        if (tempAcceptedSeeds >= (j - e)) {
            sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
            totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

            #pragma omp critical
            {
                numLocations += totalPossibleLocations.size();
                possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
            };
        } else if (tempAcceptedSeeds == 0 && isReverseAccepted) {
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
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.size() < q) {
                        startPosition += seed.size() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                exactSearchingPosition(seed, reverseRead, mode, q, k, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());
                }
            }

            if (tempAcceptedSeeds >= (j - e)) {
                sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                #pragma omp critical
                {
                    numLocations += totalPossibleLocations.size();
                    possibleReadsMap[reverseRead] = vector<unsigned long long int>(totalPossibleLocations);
                };
            }
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void approximateSearchingForAll() {
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

            approximateSearchingPosition(seed, forwardRead, mode, q, k, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());
            }
        }

        if (tempAcceptedSeeds >= (j - e % j)) {
            sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
            totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

            #pragma omp critical
            {
                numLocations += totalPossibleLocations.size();
                possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
            };
        } else if (tempAcceptedSeeds == 0 && isReverseAccepted) {
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
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.size() < q) {
                        startPosition += seed.size() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                approximateSearchingPosition(seed, reverseRead, mode, q, k, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());
                }
            }

            if (tempAcceptedSeeds >= (j - e % j)) {
                sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                #pragma omp critical
                {
                    numLocations += totalPossibleLocations.size();
                    possibleReadsMap[reverseRead] = vector<unsigned long long int>(totalPossibleLocations);
                };
            }
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void searchingReadProcess() {
    j = m / q;
    allowableE = floor(e / j);

//    bool isExactMatching = false;
//
//    if (allowableE == 0) {
//        isExactMatching = true;
//    }

    if (allowableE == 0) {
        exactSearchingForAll();
    } else {
        approximateSearchingForAll();
    }
}

void exactSearchingForExit() {
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

            exactSearchingPosition(seed, forwardRead, mode, q, k, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());

                if (tempAcceptedSeeds == (j - e)) {
                    sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                    totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                    #pragma omp critical
                    {
                        numLocations += totalPossibleLocations.size();
                        possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
                    };
                    break;
                }
            }
        }

        if (tempAcceptedSeeds == 0 && isReverseAccepted) {
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
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.size() < q) {
                        startPosition += seed.size() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                exactSearchingPosition(seed, reverseRead, mode, q, k, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());

                    if (tempAcceptedSeeds == (j - e)) {
                        sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                        totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                        #pragma omp critical
                        {
                            numLocations += totalPossibleLocations.size();
                            possibleReadsMap[reverseRead] = vector<unsigned long long int>(totalPossibleLocations);
                        };
                        break;
                    }
                }
            }
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void approximateSearchingForExit() {
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

            approximateSearchingPosition(seed, forwardRead, mode, q, k, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());

                if (tempAcceptedSeeds == (j - e % j)) {
                    sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                    totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                    #pragma omp critical
                    {
                        numLocations += totalPossibleLocations.size();
                        possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
                    };
                    break;
                }
            }
        }

        if (tempAcceptedSeeds == 0 && isReverseAccepted) {
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
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.size() < q) {
                        startPosition += seed.size() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                approximateSearchingPosition(seed, reverseRead, mode, q, k, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());

                    if (tempAcceptedSeeds == (j - e % j)) {
                        sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                        totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                        #pragma omp critical
                        {
                            numLocations += totalPossibleLocations.size();
                            possibleReadsMap[reverseRead] = vector<unsigned long long int>(totalPossibleLocations);
                        };
                        break;
                    }
                }
            }
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void searchingReadFoundExitProcess() {
    j = m / q;
    allowableE = floor(e / j);

//    bool isExactMatching = false;
//
//    if (allowableE == 0) {
//        isExactMatching = true;
//    }

    if (allowableE == 0) {
        exactSearchingForExit();
    } else {
        approximateSearchingForExit();
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
