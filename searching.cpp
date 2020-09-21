//
// Created by saimanalili on 25/02/2020.
//

#include "searching.h"
#include "seedselector.h"

vector<unsigned long long int> searchingUsingMinimizers(string seed, int adjustmentValue, vector<unsigned long long int>& location) {
    vector<unsigned long long int> foundLocations;

    int i;
    for (i = 0; i < location.size(); i++) {
        if (seed.compare(refGenome.genomeData.substr(location[i], q)) == 0) {
            if ((location[i] - adjustmentValue) >= 0 && (location[i] - adjustmentValue) < refGenome.genomeData.size()) {
                    foundLocations.push_back(location[i] - adjustmentValue);
            }
        }
    }

    return foundLocations;
}

vector<unsigned long long int> approximateSearchingUsingMinimizers(string seed, int adjustmentValue, vector<unsigned long long int>& location) {
    vector<unsigned long long int> foundLocations;
    int i;
    for (i = 0; i < location.size(); i++) {
        EdlibAlignResult result = edlibAlign(refGenome.genomeData.substr(location[i], q).c_str(), q,
                seed.c_str(), seed.length(), edlibDefaultAlignConfig());

        if (result.editDistance <= allowableE) {
            if ((location[i] - adjustmentValue) >= 0 && (location[i] - adjustmentValue) < refGenome.genomeData.size()) {
                foundLocations.push_back(location[i] - adjustmentValue);
            }
        }

        edlibFreeAlignResult(result);
    }

    return foundLocations;
}

vector<unsigned long long int> searchingUsingDirectOrOpen(string seed, unsigned long long int index, int adjustmentValue, vector<unsigned long long int>& location) {
    vector<unsigned long long int> foundLocations;

    while (seed.compare(refGenome.genomeData.substr(posTable[index], q)) == 0) {
        if ((posTable[index] - adjustmentValue) >= 0 && (posTable[index] - adjustmentValue) < refGenome.genomeData.size()) {
            foundLocations.push_back(posTable[index] - adjustmentValue);
        }

        index++;
    }

    return foundLocations;
}

vector<unsigned long long int> approximateSearchingUsingDirectOrOpen(string seed, unsigned long long int index, int adjustmentValue, vector<unsigned long long int>& location) {
    vector<unsigned long long int> foundLocations;
    bool continueCompare = true;

    while (continueCompare) {
        if (index < posTable.size()) {
            EdlibAlignResult result = edlibAlign(refGenome.genomeData.substr(posTable[index], q).c_str(),
                                                 q, seed.c_str(), seed.length(), edlibDefaultAlignConfig());

            if (result.editDistance <= allowableE) {
                if ((posTable[index] - adjustmentValue) >= 0 && (posTable[index] - adjustmentValue) < refGenome.genomeData.size()) {
                    foundLocations.push_back(posTable[index] - adjustmentValue);
                }

                index++;
            } else {
                continueCompare = false;
            }

            edlibFreeAlignResult(result);
        } else {
            continueCompare = false;
        }
    }

    return foundLocations;
}

vector<unsigned long long int> exactSearchingPosition(string seed, string mode, int adjustmentValue) {
    vector<unsigned long long int> location;
    vector<unsigned long long int> foundLocations;

    unsigned long long int rank;
    unsigned long long int index;

    if (mode.compare("min") == 0) {
        rank = getMinimizerRank(seed, q, w);
        location = minimizers[rank];

        if (location.size() > 0) {
            foundLocations = searchingUsingMinimizers(seed.substr(0, q), adjustmentValue, location);
        }
    } else if (mode.compare("dir") == 0 || mode.compare("open") == 0) {
        if (mode.compare("dir") == 0) {
            rank = extractRanking(seed);

            if (rank >= 0) {
                index = dirTable[rank];

                if (index < posTable.size()) {
                    foundLocations = searchingUsingDirectOrOpen(seed, index, adjustmentValue, location);
                }
            }
        } else {
            rank = extractRanking(seed);

            try {
                index = dirTable[codeTable[rank]];

                foundLocations = searchingUsingDirectOrOpen(seed, index, adjustmentValue, location);
            } catch (exception& e) { }
        }
    }

    return foundLocations;
}

vector<unsigned long long int> approximateSearchingPosition(string seed, string mode, int adjustmentValue) {
    vector<unsigned long long int> location;
    vector<unsigned long long int> foundLocations;

    unsigned long long int rank;
    unsigned long long int index;

    if (mode.compare("min") == 0) {
        rank = getMinimizerRank(seed, q, w);
        location = minimizers[rank];

        if (location.size() > 0) {
            foundLocations = approximateSearchingUsingMinimizers(seed.substr(0, q), adjustmentValue, location);
        }
    } else if (mode.compare("dir") == 0 || mode.compare("open") == 0) {
        if (mode.compare("dir") == 0) {
            rank = extractRanking(seed);

            if (rank >= 0) {
                index = dirTable[rank];

                if (index < posTable.size()) {
                    foundLocations = approximateSearchingUsingDirectOrOpen(seed, index, adjustmentValue, location);
                }
            }
        } else {
            rank = extractRanking(seed);

            try {
                index = dirTable[codeTable[rank]];

                foundLocations = approximateSearchingUsingDirectOrOpen(seed, index, adjustmentValue, location);
            } catch (exception& e) { }
        }
    }

    return foundLocations;
}