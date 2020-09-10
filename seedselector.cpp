//
// Created by saimanalili on 25/02/2020.
//

#include "seedselector.h"

int j;
int allowableE;

void exactSearchingForAll() {
    int acceptanceCriterion = (j - e);

    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;
            int adjustmentValue = startPosition;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);

                if (seed.length() < w) {
                    adjustmentValue = startPosition - (w - seed.length());
                    startPosition += seed.length() - w;
                    seed = forwardRead.substr(startPosition, w);
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.length() < q) {
                    adjustmentValue = startPosition - (q - seed.length());
                    startPosition += seed.length() - q;
                    seed = forwardRead.substr(startPosition, q);
                }
            }

            vector<unsigned long long int> forward;

            exactSearchingPosition(seed, mode, adjustmentValue, k, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());
            }
        }

        if (tempAcceptedSeeds >= acceptanceCriterion) {
            sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
            totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

            #pragma omp critical
            {
                numPossibleReadLocations += totalPossibleLocations.size();
                reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleLocations);
                numAcceptedReads++;
            };
        } else if (tempAcceptedSeeds < acceptanceCriterion && isReverseAccepted) {
            tempAcceptedSeeds = 0;
            string reverseRead(reverseComplement(forwardRead));

            int k;
            for (k = 0; k < j; k++) {
                int startPosition = k * q;
                int adjustmentValue = startPosition;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);

                    if (seed.length() < w) {
                        adjustmentValue = startPosition - (w - seed.length());
                        startPosition += seed.length() - w;
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.length() < q) {
                        adjustmentValue = startPosition - (q - seed.length());
                        startPosition += seed.length() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                exactSearchingPosition(seed, mode, adjustmentValue, k, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());
                }
            }

            if (tempAcceptedSeeds >= acceptanceCriterion) {
                sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                #pragma omp critical
                {
                    numPossibleReadLocations += totalPossibleLocations.size();
                    reads[i].reverseLocations = vector<unsigned long long int>(totalPossibleLocations);
                    numAcceptedReads++;
                };
            }
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void approximateSearchingForAll() {
    int acceptanceCriterion;

    if (allowableE < q) {
        acceptanceCriterion = j - floor(e / (allowableE + 1));
    } else {
        acceptanceCriterion = j - floor(e / allowableE);
    }

    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;
            int adjustmentValue = startPosition;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);

                if (seed.length() < w) {
                    adjustmentValue = startPosition - (w - seed.length());
                    startPosition += seed.length() - w;
                    seed = forwardRead.substr(startPosition, w);
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.length() < q) {
                    adjustmentValue = startPosition - (q - seed.length());
                    startPosition += seed.length() - q;
                    seed = forwardRead.substr(startPosition, q);
                }
            }

            vector<unsigned long long int> forward;

            approximateSearchingPosition(seed, mode, adjustmentValue, k, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());
            }
        }

        if (tempAcceptedSeeds >= acceptanceCriterion) {
            sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
            totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

            #pragma omp critical
            {
                numPossibleReadLocations += totalPossibleLocations.size();
                reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleLocations);
                numAcceptedReads++;
            };
        } else if (tempAcceptedSeeds < acceptanceCriterion && isReverseAccepted) {
            tempAcceptedSeeds = 0;
            string reverseRead(reverseComplement(forwardRead));

            int k;
            for (k = 0; k < j; k++) {
                int startPosition = k * q;
                int adjustmentValue = startPosition;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);

                    if (seed.length() < w) {
                        adjustmentValue = startPosition - (w - seed.length());
                        startPosition += seed.length() - w;
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.length() < q) {
                        adjustmentValue = startPosition - (q - seed.length());
                        startPosition += seed.length() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                approximateSearchingPosition(seed, mode, adjustmentValue, k, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());
                }
            }

            if (tempAcceptedSeeds >= acceptanceCriterion) {
                sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                #pragma omp critical
                {
                    numPossibleReadLocations += totalPossibleLocations.size();
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
    int acceptanceCriterion = (j - e);

    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;
            int adjustmentValue = startPosition;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);

                if (seed.length() < w) {
                    adjustmentValue = startPosition - (w - seed.length());
                    startPosition += seed.length() - w;
                    seed = forwardRead.substr(startPosition, w);
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.length() < q) {
                    adjustmentValue = startPosition - (q - seed.length());
                    startPosition += seed.length() - q;
                    seed = forwardRead.substr(startPosition, q);
                }
            }

            vector<unsigned long long int> forward;

            exactSearchingPosition(seed, mode, adjustmentValue, k, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());

                if (tempAcceptedSeeds == acceptanceCriterion) {
                    sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                    totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                    #pragma omp critical
                    {
                        numPossibleReadLocations += totalPossibleLocations.size();
                        reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleLocations);
                        numAcceptedReads++;
                    };
                    break;
                }
            }
        }

        if (tempAcceptedSeeds < acceptanceCriterion && isReverseAccepted) {
            tempAcceptedSeeds = 0;
            string reverseRead(reverseComplement(forwardRead));

            int k;
            for (k = 0; k < j; k++) {
                int startPosition = k * q;
                int adjustmentValue = startPosition;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);

                    if (seed.length() < w) {
                        adjustmentValue = startPosition - (w - seed.length());
                        startPosition += seed.length() - w;
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.length() < q) {
                        adjustmentValue = startPosition - (q - seed.length());
                        startPosition += seed.length() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                exactSearchingPosition(seed, mode, adjustmentValue, k, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());

                    if (tempAcceptedSeeds == acceptanceCriterion) {
                        sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                        totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                        #pragma omp critical
                        {
                            numPossibleReadLocations += totalPossibleLocations.size();
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
    int acceptanceCriterion;

    if (allowableE < q) {
        acceptanceCriterion = j - floor(e / (allowableE + 1));
    } else {
        acceptanceCriterion = j - floor(e / allowableE);
    }

    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;
            int adjustmentValue = startPosition;

            string seed;

            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, w);

                if (seed.length() < w) {
                    adjustmentValue = startPosition - (w - seed.length());
                    startPosition += seed.length() - w;
                    seed = forwardRead.substr(startPosition, w);
                }
            } else {
                seed = forwardRead.substr(startPosition, q);

                if (seed.length() < q) {
                    adjustmentValue = startPosition - (q - seed.length());
                    startPosition += seed.length() - q;
                    seed = forwardRead.substr(startPosition, q);
                }
            }

            vector<unsigned long long int> forward;

            approximateSearchingPosition(seed, mode, adjustmentValue, k, forward);

            if (forward.size() > 0) {
                tempAcceptedSeeds++;
                totalPossibleLocations.insert(totalPossibleLocations.end(), forward.begin(), forward.end());

                if (tempAcceptedSeeds >= acceptanceCriterion) {
                    sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                    totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                    #pragma omp critical
                    {
                        numPossibleReadLocations += totalPossibleLocations.size();
                        reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleLocations);
                    };
                    break;
                }
            }
        }

        if (tempAcceptedSeeds < acceptanceCriterion && isReverseAccepted) {
            tempAcceptedSeeds = 0;
            string reverseRead(reverseComplement(forwardRead));

            int k;
            for (k = 0; k < j; k++) {
                int startPosition = k * q;
                int adjustmentValue = startPosition;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, w);

                    if (seed.length() < w) {
                        adjustmentValue = startPosition - (w - seed.length());
                        startPosition += seed.length() - w;
                        seed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    seed = reverseRead.substr(startPosition, q);

                    if (seed.length() < q) {
                        adjustmentValue = startPosition - (q - seed.length());
                        startPosition += seed.length() - q;
                        seed = reverseRead.substr(startPosition, q);
                    }
                }

                vector<unsigned long long int> reverse;

                approximateSearchingPosition(seed, mode, adjustmentValue, k, reverse);

                if (reverse.size() > 0) {
                    tempAcceptedSeeds++;
                    totalPossibleLocations.insert(totalPossibleLocations.end(), reverse.begin(), reverse.end());

                    if (tempAcceptedSeeds >= acceptanceCriterion) {
                        sort(totalPossibleLocations.begin(), totalPossibleLocations.end());
                        totalPossibleLocations.erase(unique(totalPossibleLocations.begin(), totalPossibleLocations.end()), totalPossibleLocations.end());

                        #pragma omp critical
                        {
                            numPossibleReadLocations += totalPossibleLocations.size();
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