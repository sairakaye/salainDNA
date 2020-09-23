/**
 * @file searchall.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the search all mode implementation of the Seed Selector.
 */

#include "searchall.h"

void exactSearchingForAll() {
    int acceptanceCriterion = (j - e);

    /*
    exactForwardSearchingAll(acceptanceCriterion);

    if (isReverseAccepted) {
        exactReverseSearchingAll(acceptanceCriterion);
    }
    */

    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead = reads[i].readData;
        string reverseRead;

        if (isReverseAccepted) {
            reverseRead = reverseComplement(reads[i].readData);
        }

        vector<unsigned long long int> totalPossibleFwdLocations;
        vector<unsigned long long int> totalPossibleRevLocations;
        unsigned int tempAcceptedFwdSeeds = 0;
        unsigned int tempAcceptedRevSeeds = 0;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;
            int adjustmentValue = startPosition;

            string fwdSeed;
            string revSeed;

            if (mode.compare("min") == 0) {
                fwdSeed = forwardRead.substr(startPosition, w);

                if (fwdSeed.length() < w) {
                    adjustmentValue = startPosition - (w - fwdSeed.length());
                    startPosition += fwdSeed.length() - w;
                    fwdSeed = forwardRead.substr(startPosition, w);
                }
            } else {
                fwdSeed = forwardRead.substr(startPosition, q);

                if (fwdSeed.length() < q) {
                    adjustmentValue = startPosition - (q - fwdSeed.length());
                    startPosition += fwdSeed.length() - q;
                    fwdSeed = forwardRead.substr(startPosition, q);
                }
            }

            if (isReverseAccepted) {
                if (mode.compare("min") == 0) {
                    revSeed = reverseRead.substr(startPosition, w);

                    if (revSeed.length() < w) {
                        adjustmentValue = startPosition - (w - revSeed.length());
                        startPosition += revSeed.length() - w;
                        revSeed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    revSeed = reverseRead.substr(startPosition, q);

                    if (revSeed.length() < q) {
                        adjustmentValue = startPosition - (q - revSeed.length());
                        startPosition += revSeed.length() - q;
                        revSeed = reverseRead.substr(startPosition, q);
                    }
                }
            }

            vector<unsigned long long int> forward;

            forward = exactSearchingPosition(fwdSeed, mode, adjustmentValue);

            if (forward.size() > 0) {
                tempAcceptedFwdSeeds++;
                totalPossibleFwdLocations.insert(totalPossibleFwdLocations.end(), forward.begin(), forward.end());
            }

            if (isReverseAccepted) {
                vector<unsigned long long int> reverse;

                reverse = exactSearchingPosition(revSeed, mode, adjustmentValue);

                if (reverse.size() > 0) {
                    tempAcceptedRevSeeds++;
                    totalPossibleRevLocations.insert(totalPossibleRevLocations.end(), reverse.begin(), reverse.end());
                }
            }
        }

        if (tempAcceptedFwdSeeds >= acceptanceCriterion) {
            sort(totalPossibleFwdLocations.begin(), totalPossibleFwdLocations.end());
            totalPossibleFwdLocations.erase(unique(totalPossibleFwdLocations.begin(), totalPossibleFwdLocations.end()), totalPossibleFwdLocations.end());

            #pragma omp critical
            {
                numPossibleReadLocations += totalPossibleFwdLocations.size();
                reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleFwdLocations);
                numAcceptedReads++;
            };
        }

        if (isReverseAccepted) {
            if (tempAcceptedRevSeeds >= acceptanceCriterion) {
                sort(totalPossibleRevLocations.begin(), totalPossibleRevLocations.end());
                totalPossibleRevLocations.erase(unique(totalPossibleRevLocations.begin(), totalPossibleRevLocations.end()), totalPossibleRevLocations.end());

                #pragma omp critical
                {
                    numPossibleReadLocations += totalPossibleRevLocations.size();
                    reads[i].reverseLocations = vector<unsigned long long int>(totalPossibleRevLocations);
                    numAcceptedReads++;
                };
            }
        }

        numAcceptedSeeds += tempAcceptedFwdSeeds + tempAcceptedRevSeeds;
    }
}

void approximateSearchingForAll() {
    int acceptanceCriterion;

    if (allowableE < q) {
        acceptanceCriterion = j - floor(e / (allowableE + 1));
    } else {
        acceptanceCriterion = j - floor(e / allowableE);
    }

    /*
    approximateForwardSearchingAll(acceptanceCriterion);

    if (isReverseAccepted) {
        approximateReverseSearchingAll(acceptanceCriterion);
    }
    */

    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead = reads[i].readData;
        string reverseRead;

        if (isReverseAccepted) {
            reverseRead = reverseComplement(reads[i].readData);
        }

        vector<unsigned long long int> totalPossibleFwdLocations;
        vector<unsigned long long int> totalPossibleRevLocations;
        unsigned int tempAcceptedFwdSeeds = 0;
        unsigned int tempAcceptedRevSeeds = 0;

        int k;
        for (k = 0; k < j; k++) {
            int startPosition = k * q;
            int adjustmentValue = startPosition;

            string fwdSeed;
            string revSeed;

            if (mode.compare("min") == 0) {
                fwdSeed = forwardRead.substr(startPosition, w);

                if (fwdSeed.length() < w) {
                    adjustmentValue = startPosition - (w - fwdSeed.length());
                    startPosition += fwdSeed.length() - w;
                    fwdSeed = forwardRead.substr(startPosition, w);
                }
            } else {
                fwdSeed = forwardRead.substr(startPosition, q);

                if (fwdSeed.length() < q) {
                    adjustmentValue = startPosition - (q - fwdSeed.length());
                    startPosition += fwdSeed.length() - q;
                    fwdSeed = forwardRead.substr(startPosition, q);
                }
            }

            if (isReverseAccepted) {
                if (mode.compare("min") == 0) {
                    revSeed = reverseRead.substr(startPosition, w);

                    if (revSeed.length() < w) {
                        adjustmentValue = startPosition - (w - revSeed.length());
                        startPosition += revSeed.length() - w;
                        revSeed = reverseRead.substr(startPosition, w);
                    }
                } else {
                    revSeed = reverseRead.substr(startPosition, q);

                    if (revSeed.length() < q) {
                        adjustmentValue = startPosition - (q - revSeed.length());
                        startPosition += revSeed.length() - q;
                        revSeed = reverseRead.substr(startPosition, q);
                    }
                }
            }

            vector<unsigned long long int> forward;

            forward = approximateSearchingPosition(fwdSeed, mode, adjustmentValue);

            if (forward.size() > 0) {
                tempAcceptedFwdSeeds++;
                totalPossibleFwdLocations.insert(totalPossibleFwdLocations.end(), forward.begin(), forward.end());
            }

            if (isReverseAccepted) {
                vector<unsigned long long int> reverse;

                reverse = approximateSearchingPosition(revSeed, mode, adjustmentValue);

                if (reverse.size() > 0) {
                    tempAcceptedRevSeeds++;
                    totalPossibleRevLocations.insert(totalPossibleRevLocations.end(), reverse.begin(), reverse.end());
                }
            }
        }

        if (tempAcceptedFwdSeeds >= acceptanceCriterion) {
            sort(totalPossibleFwdLocations.begin(), totalPossibleFwdLocations.end());
            totalPossibleFwdLocations.erase(unique(totalPossibleFwdLocations.begin(), totalPossibleFwdLocations.end()), totalPossibleFwdLocations.end());

            #pragma omp critical
            {
                numPossibleReadLocations += totalPossibleFwdLocations.size();
                reads[i].forwardLocations = vector<unsigned long long int>(totalPossibleFwdLocations);
                numAcceptedReads++;
            };
        }

        if (isReverseAccepted) {
            if (tempAcceptedRevSeeds >= acceptanceCriterion) {
                sort(totalPossibleRevLocations.begin(), totalPossibleRevLocations.end());
                totalPossibleRevLocations.erase(unique(totalPossibleRevLocations.begin(), totalPossibleRevLocations.end()), totalPossibleRevLocations.end());

                #pragma omp critical
                {
                    numPossibleReadLocations += totalPossibleRevLocations.size();
                    reads[i].reverseLocations = vector<unsigned long long int>(totalPossibleRevLocations);
                    numAcceptedReads++;
                };
            }
        }

        numAcceptedSeeds += tempAcceptedFwdSeeds + tempAcceptedRevSeeds;
    }
}

/*
void exactForwardSearchingAll(int acceptanceCriterion) {
    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead = reads[i].readData;
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
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void exactReverseSearchingAll(int acceptanceCriterion) {
    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string reverseRead = reverseComplement(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

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

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}
*/

/*
void approximateForwardSearchingAll(int acceptanceCriterion) {
    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string forwardRead = reads[i].readData;
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
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void approximateReverseSearchingAll(int acceptanceCriterion) {
    int i;
    #pragma omp parallel for reduction(+:numAcceptedSeeds)
    for (i = 0; i < reads.size(); i++) {
        string reverseRead = reverseComplement(reads[i].readData);
        vector<unsigned long long int> totalPossibleLocations;
        unsigned int tempAcceptedSeeds = 0;

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

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}
*/