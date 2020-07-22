//
// Created by saimanalili on 25/02/2020.
//

#include "pigeonhole.h"

int j;
int allowableE;

void exactSearchingForAll() {
    /* For DEBUGGING only.
    ofstream infoFile;
    string infoFileName("testing_hatdog.txt");
    infoFile.open(infoFileName.c_str(), ios::out);
    */

    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);

                if (seed.length() < w) {
                    startPosition += seed.length() - w;
                    seed = forwardRead.substr(startPosition, w);
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.length() < q) {
                    startPosition += seed.length() - q;
                    seed = forwardRead.substr(startPosition, q);
                }
            }

            vector<unsigned long long int> forward;

            exactSearchingPosition(seed, mode, q, k, forward);

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
                numPossibleReadLocations += totalPossibleLocations.size();
                //possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
                reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleLocations);
                numAcceptedReads++;
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

                    if (seed.length() < w) {
                        startPosition += seed.length() - w;
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.length() < q) {
                        startPosition += seed.length() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                exactSearchingPosition(seed, mode, q, k, reverse);

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
                    numPossibleReadLocations += totalPossibleLocations.size();
                    //possibleReadsMap[reverseRead] = vector<unsigned long long int>(totalPossibleLocations);
                    reads[i].reverseLocations = vector<unsigned long long int>(totalPossibleLocations);
                    numAcceptedReads++;
                };
            }
        } else {
            /*
            #pragma omp critical
            infoFile << ">" << reads[i].readName << "\n" << reads[i].readData << endl;
            */
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }

    //infoFile.close();
}

void approximateSearchingForAll() {
    int i;
    //#pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

        int minEditFound = INT_MAX;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);

                if (seed.length() < w) {
                    startPosition += seed.length() - w;
                    seed = forwardRead.substr(startPosition, w);
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.length() < q) {
                    startPosition += seed.length() - q;
                    seed = forwardRead.substr(startPosition, q);
                }
            }

            vector<unsigned long long int> forward;

            approximateSearchingPosition(seed, mode, q, k, forward, &minEditFound);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());
            }
        }

        if (tempAcceptedSeeds == j && minEditFound <= allowableE) {
            sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
            totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

            #pragma omp critical
            {
                numPossibleReadLocations += totalPossibleLocations.size();
                //possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
                reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleLocations);
                numAcceptedReads++;
            };
        } else if (tempAcceptedSeeds == 0 && isReverseAccepted) {
            tempAcceptedSeeds = 0;
            string reverseRead = reverseComplement(forwardRead);

            int minEditFound = INT_MAX;

            int k;
            for (k = 0; k < j; k++) {
                int startPosition = k * q;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);

                    if (seed.length() < w) {
                        startPosition += seed.length() - w;
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.length() < q) {
                        startPosition += seed.length() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                approximateSearchingPosition(seed, mode, q, k, reverse, &minEditFound);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());
                }
            }

            if (tempAcceptedSeeds == j && minEditFound <= allowableE) {
                sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                #pragma omp critical
                {
                    numPossibleReadLocations += totalPossibleLocations.size();
                    //possibleReadsMap[reverseRead] = vector<unsigned long long int>(totalPossibleLocations);
                    reads[i].reverseLocations = vector<unsigned long long int>(totalPossibleLocations);
                    numAcceptedReads++;
                };
            }
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void searchingReadProcess() {
    j = ceil(m / (double) q);
    allowableE = floor(e / (double) j);

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
        string forwardRead(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);

                if (seed.length() < w) {
                    startPosition += seed.length() - w;
                    seed = forwardRead.substr(startPosition, w);
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.length() < q) {
                    startPosition += seed.length() - q;
                    seed = forwardRead.substr(startPosition, q);
                }
            }

            vector<unsigned long long int> forward;

            exactSearchingPosition(seed, mode, q, k, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());

                if (tempAcceptedSeeds == (j - e)) {
                    sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                    totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                    #pragma omp critical
                    {
                        numPossibleReadLocations += totalPossibleLocations.size();
                        //possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
                        reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleLocations);
                        numAcceptedReads++;
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

                    if (seed.length() < w) {
                        startPosition += seed.length() - w;
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.length() < q) {
                        startPosition += seed.length() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                exactSearchingPosition(seed, mode, q, k, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());

                    if (tempAcceptedSeeds == (j - e)) {
                        sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                        totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                        #pragma omp critical
                        {
                            numPossibleReadLocations += totalPossibleLocations.size();
                            //possibleReadsMap[reverseRead] = vector<unsigned long long int>(totalPossibleLocations);
                            reads[i].reverseLocations = vector<unsigned long long int>(totalPossibleLocations);
                            numAcceptedReads++;
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
        string forwardRead(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

        int minEditFound = INT_MAX;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);

                if (seed.length() < w) {
                    startPosition += seed.length() - w;
                    seed = forwardRead.substr(startPosition, w);
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.length() < q) {
                    startPosition += seed.length() - q;
                    seed = forwardRead.substr(startPosition, q);
                }
            }

            vector<unsigned long long int> forward;

            approximateSearchingPosition(seed, mode, q, k, forward, &minEditFound);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());

                if (tempAcceptedSeeds == (j - e % j)) {
                    sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                    totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                    #pragma omp critical
                    {
                        numPossibleReadLocations += totalPossibleLocations.size();
                        //possibleReadsMap[forwardRead] = vector<unsigned long long int>(totalPossibleLocations);
                        reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleLocations);
                    };
                    break;
                }
            }
        }

        if (tempAcceptedSeeds == 0 && isReverseAccepted) {
            tempAcceptedSeeds = 0;
            string reverseRead = reverseComplement(forwardRead);

            int minEditFound = INT_MAX;

            int k;
            for (k = 0; k < j; k++) {
                int startPosition = k * q;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);

                    if (seed.length() < w) {
                        startPosition += seed.length() - w;
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.length() < q) {
                        startPosition += seed.length() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                approximateSearchingPosition(seed, mode, q, k, reverse, &minEditFound);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());

                    if (tempAcceptedSeeds == (j - e % j)) {
                        sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                        totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                        #pragma omp critical
                        {
                            numPossibleReadLocations += totalPossibleLocations.size();
                            //possibleReadsMap[reverseRead] = vector<unsigned long long int>(totalPossibleLocations);
                            reads[i].reverseLocations = vector<unsigned long long int>(totalPossibleLocations);
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
    j = ceil(m / (double) q);
    allowableE = floor(e / (double) j);

    if (allowableE == 0) {
        exactSearchingForExit();
    } else {
        approximateSearchingForExit();
    }
}