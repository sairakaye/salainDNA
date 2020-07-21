//
// Created by saimanalili on 25/02/2020.
//

#include "searching.h"
#include "pigeonhole.h"

void searchingUsingMinimizers(string seed, string read, int k, vector<unsigned long long int>& foundLocations,
vector<unsigned long long int>& location) {
    int i;
    for (i = 0; i < location.size(); i++) {
        if (seed.compare(refGenome.genomeData.substr(location[i], q)) == 0) {
            if ((location[i] - (q * k)) >= 0 && (location[i] - (q * k)) < refGenome.genomeData.size()) {
                //#pragma omp critical
                //{
                    foundLocations.push_back(location[i] - (q * k));

                    //if (isForwardStrand) {
                    //    forwardReadsMap[read].push_back(location[i] - (q * k));
                    //} else {
                    //    reverseReadsMap[read].push_back(location[i] - (q * k));
                   //}
                //};
            }
        }
    }
}

void approximateSearchingUsingMinimizers(string seed, string read, int k, vector<unsigned long long int>& foundLocations,
vector<unsigned long long int>& location) {
    int i;
    for (i = 0; i < location.size(); i++) {
        EdlibAlignResult result = edlibAlign(refGenome.genomeData.substr(location[i], seed.length()).c_str(), seed.length(), seed.c_str(), seed.length(), edlibDefaultAlignConfig());

        if (result.editDistance <= allowableE) {
            if ((location[i] - (q * k)) >= 0 && (location[i] - (q * k)) < refGenome.genomeData.size()) {
                //#pragma omp critical
                //{
                foundLocations.push_back(location[i] - (q * k));

                //if (isForwardStrand) {
                //    forwardReadsMap[read].push_back(location[i] - (q * k));
                //} else {
                //    reverseReadsMap[read].push_back(location[i] - (q * k));
                //}
                //};
            }
        }

        edlibFreeAlignResult(result);
    }
}

void searchingUsingDirectOrOpen(string seed, string read, unsigned long long int index, string mode, int k, vector<unsigned long long int>& foundLocations,
vector<unsigned long long int>& location) {
    while (seed.compare(refGenome.genomeData.substr(posTable[index], seed.size())) == 0) {
        if ((posTable[index] - (q * k)) >= 0 && (posTable[index]- (q * k)) < refGenome.genomeData.size()) {
            //#pragma omp critical
            //{
            foundLocations.push_back(posTable[index] - (q * k));

            //if (isForwardStrand) {
            //    forwardReadsMap[read].push_back(posTable[index]  - (q * k));
            //} else {
            //    reverseReadsMap[read].push_back(posTable[index]  - (q * k));
            //}
            //};
        }

        index++;
    }
}

void approximateSearchingUsingDirectOrOpen(string seed, string read, unsigned long long int index, string mode, int k, vector<unsigned long long int>& foundLocations,
vector<unsigned long long int>& location) {
    bool continueCompare = true;

    while (continueCompare) {
        EdlibAlignResult result = edlibAlign(refGenome.genomeData.substr(posTable[index], seed.length()).c_str(), seed.length(), seed.c_str(), seed.length(), edlibDefaultAlignConfig());

        if (result.editDistance <= allowableE) {
            if ((posTable[index] - (q * k)) >= 0 && (posTable[index]- (q * k)) < refGenome.genomeData.size()) {
                //#pragma omp critical
                //{
                foundLocations.push_back(posTable[index] - (q * k));

                //if (isForwardStrand) {
                //    forwardReadsMap[read].push_back(posTable[index]  - (q * k));
                //} else {
                //    reverseReadsMap[read].push_back(posTable[index]  - (q * k));
                //}
                //};
            }

            index++;
        } else {
            continueCompare = false;
        }

        edlibFreeAlignResult(result);
    }
}

void exactSearchingPosition(string seed, string read, string mode, int q, int k, vector<unsigned long long int>& foundLocations) {
    vector<unsigned long long int> location;

    unsigned long long int rank;
    unsigned long long int index;

    if (mode.compare("min") == 0) {
        rank = getMinimizerRank(seed, q, w);
        location = minimizers[rank];

        if (location.size() > 0) {
            searchingUsingMinimizers(seed.substr(0, q), read, k, foundLocations, location);
        }
    } else if (mode.compare("dir") == 0 || mode.compare("open") == 0) {
        if (mode.compare("dir") == 0) {
            rank = extractRanking(seed);

            if (rank >= 0) {
                index = dirTable[rank];

                if (index < posTable.size()) {
                    searchingUsingDirectOrOpen(seed, read, index, mode, k, foundLocations, location);
                }
            }
        } else {
            rank = extractRanking(seed);

            try {
                index = dirTable[codeTable[rank]];

                searchingUsingDirectOrOpen(seed, read, index, mode, k, foundLocations, location);
            } catch (exception& e) { }
        }
    }
}

void approximateSearchingPosition(string seed, string read, string mode, int q, int k, vector<unsigned long long int>& foundLocations) {
    vector<unsigned long long int> location;

    unsigned long long int rank;
    unsigned long long int index;

    if (mode.compare("min") == 0) {
        rank = getMinimizerRank(seed, q, w);
        location = minimizers[rank];

        if (location.size() > 0) {
            approximateSearchingUsingMinimizers(seed.substr(0, q), read, k, foundLocations, location);
        }
    } else if (mode.compare("dir") == 0 || mode.compare("open") == 0) {
        if (mode.compare("dir") == 0) {
            rank = extractRanking(seed);

            if (rank >= 0) {
                index = dirTable[rank];

                if (index < posTable.size()) {
                    approximateSearchingUsingDirectOrOpen(seed, read, index, mode, k, foundLocations, location);
                }
            }
        } else {
            rank = extractRanking(seed);

            try {
                index = dirTable[codeTable[rank]];

                approximateSearchingUsingDirectOrOpen(seed, read, index, mode, k,foundLocations, location);
            } catch (exception& e) { }
        }
    }
}