/**
 * @file searchexit.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the exit mode implementation of the Seed Selector.
 */

#include "searchexit.h"

void exactForwardSearchingExit(int acceptanceCriterion) {
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

            forward = exactSearchingPosition(seed, mode, adjustmentValue);

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

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void exactReverseSearchingExit(int acceptanceCriterion) {
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

            reverse = exactSearchingPosition(seed, mode, adjustmentValue);

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

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void exactSearchingForExit() {
    int acceptanceCriterion = (j - e);

    exactForwardSearchingExit(acceptanceCriterion);

    if (isReverseAccepted) {
        exactReverseSearchingExit(acceptanceCriterion);
    }
}

void approximateForwardSearchingExit(int acceptanceCriterion) {
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

            forward = approximateSearchingPosition(seed, mode, adjustmentValue);

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
                        numAcceptedReads++;
                    };

                    break;
                }
            }
        }

        numAcceptedSeeds += tempAcceptedSeeds;
    }
}

void approximateReverseSearchingExit(int acceptanceCriterion) {
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

            reverse = approximateSearchingPosition(seed, mode, adjustmentValue);

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
                        numAcceptedReads++;
                    };
                    break;
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

    approximateForwardSearchingExit(acceptanceCriterion);

    if (isReverseAccepted) {
        approximateReverseSearchingExit(acceptanceCriterion);
    }
}