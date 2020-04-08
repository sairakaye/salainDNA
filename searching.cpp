//
// Created by saimanalili on 25/02/2020.
//

#include "searching.h"
#include "pigeonhole.h"

void searchingUsingMinimizers(string seed, string read, int k, bool isForwardStrand, vector<unsigned long long int>& foundLocations,
vector<unsigned long long int>& location) {
    for (int i = 0; i < location.size(); i++) {
        if (seed.compare(refGenome.substr(location[i], q)) == 0) {
            if ((location[i] - (q * k)) >= 0 && (location[i] - (q * k)) < refGenome.size()) {
                #pragma omp critical
                {
                    foundLocations.push_back(location[i] - (q * k));

                    if (isForwardStrand) {
                        forwardReadsMap[read].push_back(location[i] - (q * k));
                    } else {
                        reverseReadsMap[read].push_back(location[i] - (q * k));
                    }
                };
            } else {
                continue;
            }
        }
    }
}

void approximateSearchingUsingMinimizers(string seed, string read, int k, bool isForwardStrand, vector<unsigned long long int>& foundLocations,
vector<unsigned long long int>& location) {
    for (int i = 0; i < location.size(); i++) {
        EdlibAlignResult result = edlibAlign(refGenome.substr(location[i], seed.length()).c_str(), seed.length(), seed.c_str(), seed.length(), edlibDefaultAlignConfig());

        if (result.status == EDLIB_STATUS_OK) {
            if (result.editDistance <= allowableE) {
                if ((location[i] - (q * k)) >= 0 && (location[i] - (q * k)) < refGenome.size()) {
                    #pragma omp critical
                    {
                    foundLocations.push_back(location[i] - (q * k));

                    if (isForwardStrand) {
                        forwardReadsMap[read].push_back(location[i] - (q * k));
                    } else {
                        reverseReadsMap[read].push_back(location[i] - (q * k));
                    }
                    };
                } else {
                    continue;
                }
            }

            edlibFreeAlignResult(result);
        }
    }
}

void searchingUsingDirectOrOpen(string seed, string read, unsigned long long int index, string mode, int k, bool isForwardStrand, vector<unsigned long long int>& foundLocations,
vector<unsigned long long int>& location) {
    while (seed.compare(refGenome.substr(posTable[index], q)) == 0) {
        if ((posTable[index] - (q * k)) >= 0 && (posTable[index]- (q * k)) < refGenome.size()) {
            #pragma omp critical
            {
            foundLocations.push_back(posTable[index] - (q * k));

            if (isForwardStrand) {
                forwardReadsMap[read].push_back(posTable[index]  - (q * k));
            } else {
                reverseReadsMap[read].push_back(posTable[index]  - (q * k));
            }
            };
        }

        index++;
    }
}

void approximateSearchingUsingDirectOrOpen(string seed, string read, unsigned long long int index, string mode, int k, bool isForwardStrand, vector<unsigned long long int>& foundLocations,
vector<unsigned long long int>& location) {
    bool continueCompare = true;

    while (continueCompare) {
        EdlibAlignResult result = edlibAlign(refGenome.substr(posTable[index], seed.length()).c_str(), seed.length(), seed.c_str(), seed.length(), edlibDefaultAlignConfig());

        if (result.status == EDLIB_STATUS_OK) {
            if (result.editDistance <= allowableE) {
                if ((posTable[index] - (q * k)) >= 0 && (posTable[index]- (q * k)) < refGenome.size()) {
                    #pragma omp critical
                    {
                    foundLocations.push_back(posTable[index] - (q * k));

                    if (isForwardStrand) {
                        forwardReadsMap[read].push_back(posTable[index]  - (q * k));
                    } else {
                        reverseReadsMap[read].push_back(posTable[index]  - (q * k));
                    }
                    };
                }

                index++;
            } else {
                continueCompare = false;
            }

            edlibFreeAlignResult(result);
        }
    }
}

void searchingPosition(string seed, string read, string mode, int q, int k, bool isForwardStrand, bool isCheckApproximate, vector<unsigned long long int>& foundLocations) {
    vector<unsigned long long int> location;

    unsigned long long int rank;
    unsigned long long int index;

    if (mode.compare("min") == 0 && seed.size() == w) {
        rank = getMinimizerRank(seed, q, w);
        location = minimizers[rank];

        if (location.size() > 0) {
            if (isCheckApproximate) {
                approximateSearchingUsingMinimizers(seed.substr(0, q), read, k, isForwardStrand, foundLocations, location);
            } else {
                searchingUsingMinimizers(seed.substr(0, q), read, k, isForwardStrand, foundLocations, location);
            }
        }
    } else if ((mode.compare("dir") == 0 || mode.compare("open") == 0) && seed.size() == q) {
        if (mode.compare("dir") == 0) {
            rank = extractRanking(seed);

            if (rank >= 0) {
                index = dirTable[rank];

                if (index < posTable.size()) {
                    if (isCheckApproximate) {
                        approximateSearchingUsingDirectOrOpen(seed, read, index, mode, k, isForwardStrand, foundLocations, location);
                    } else {
                        searchingUsingDirectOrOpen(seed, read, index, mode, k, isForwardStrand, foundLocations, location);
                    }

                }
            }
        } else {
            rank = extractRanking(seed);

            try {
                index = dirTable[codeTable[rank]];

                if (isCheckApproximate) {
                    approximateSearchingUsingDirectOrOpen(seed, read, index, mode, k, isForwardStrand, foundLocations, location);
                } else {
                    searchingUsingDirectOrOpen(seed, read, index, mode, k, isForwardStrand, foundLocations, location);
                }


            } catch (exception& e) { }
        }
    }
}