/**
 * @file searching.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the location searching implementation of the Seed Selector.
 */

#include "searching.h"
#include "seedselector.h"

/**
* It is used for minimizers to search for the locations of the seed using exact searching.
*
* @param seed - The seed that will be used for determining its locations in the reference genome.
* @param adjustmentValue - The value for the alignment of the locations found of the seed to the read.
* @param minRankLocations - A list of locations found from the minimizer rank.
* @return found locations of the seed
*/
vector<unsigned long long int> exactSearchingUsingMinimizers(string seed, int adjustmentValue, vector<unsigned long long int>& minRankLocations) {
    vector<unsigned long long int> foundLocations;

    int i;
    for (i = 0; i < minRankLocations.size(); i++) {
        if (seed.compare(refGenome.genomeData.substr(minRankLocations[i], q)) == 0) {
            if ((minRankLocations[i] - adjustmentValue) >= 0 && (minRankLocations[i] - adjustmentValue) < refGenome.genomeData.size()) {
                    foundLocations.push_back(minRankLocations[i] - adjustmentValue);
            }
        }
    }

    return foundLocations;
}

/**
* It is used for minimizers to search for the locations of the seed using approximate searching.
*
* @param seed - The seed that will be used for determining its locations in the reference genome.
* @param adjustmentValue - The value for the alignment of the locations found of the seed to the read.
* @param minRankLocations - A list of locations found from the minimizer rank.
* @return found locations of the seed
*/
vector<unsigned long long int> approximateSearchingUsingMinimizers(string seed, int adjustmentValue, vector<unsigned long long int>& minRankLocations) {
    vector<unsigned long long int> foundLocations;
    int i;
    for (i = 0; i < minRankLocations.size(); i++) {
        EdlibAlignResult result = edlibAlign(refGenome.genomeData.substr(minRankLocations[i], q).c_str(), q,
                                             seed.c_str(), seed.length(), edlibDefaultAlignConfig());

        if (result.editDistance <= allowableE) {
            if ((minRankLocations[i] - adjustmentValue) >= 0 && (minRankLocations[i] - adjustmentValue) < refGenome.genomeData.size()) {
                foundLocations.push_back(minRankLocations[i] - adjustmentValue);
            }
        }

        edlibFreeAlignResult(result);
    }

    return foundLocations;
}

/**
* It is used for direct addressing or open addressing to search for the locations of the seed using exact searching.
*
* @param seed - The seed that will be used for determining its locations in the reference genome.
* @param index - The starting index from the position table to start the searching.
* @param adjustmentValue - The value for the alignment of the locations found of the seed to the read.
* @return found locations of the seed
*/
vector<unsigned long long int> exactSearchingUsingDirectOrOpen(string seed, unsigned long long int index, int adjustmentValue) {
    vector<unsigned long long int> foundLocations;

    while (seed.compare(refGenome.genomeData.substr(posTable[index], q)) == 0) {
        if ((posTable[index] - adjustmentValue) >= 0 && (posTable[index] - adjustmentValue) < refGenome.genomeData.size()) {
            foundLocations.push_back(posTable[index] - adjustmentValue);
        }

        index++;
    }

    return foundLocations;
}

/**
* It is used for direct addressing or open addressing to search for the locations of the seed using approximate searching.
*
* @param seed - The seed that will be used for determining its locations in the reference genome.
* @param index - The starting index from the position table to start the searching.
* @param adjustmentValue - The value for the alignment of the locations found of the seed to the read.
* @return found locations of the seed
*/
vector<unsigned long long int> approximateSearchingUsingDirectOrOpen(string seed, unsigned long long int index, int adjustmentValue) {
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

/**
* It gives the locations of the seed using exact searching.
*
* @param seed - The seed that will be used for determining its locations in the reference genome.
* @param mode - It indicates what mode of hash-based indexer will be used.
* @param adjustmentValue - The value for the alignment of the locations found of the seed to the read.
* @return found locations of the seed
*/
vector<unsigned long long int> exactSearchingPosition(string seed, string mode, int adjustmentValue) {
    vector<unsigned long long int> foundLocations;

    unsigned long long int rank;
    unsigned long long int index;

    if (mode.compare("min") == 0) {
        rank = getMinimizerRank(seed, q, w);

        vector<unsigned long long int> minRankLocations;
        minRankLocations = minimizers[rank];

        if (minRankLocations.size() > 0) {
            foundLocations = exactSearchingUsingMinimizers(seed.substr(0, q), adjustmentValue, minRankLocations);
        }
    } else if (mode.compare("dir") == 0 || mode.compare("open") == 0) {
        if (mode.compare("dir") == 0) {
            rank = extractRanking(seed);

            if (rank >= 0) {
                index = dirTable[rank];

                if (index < posTable.size()) {
                    foundLocations = exactSearchingUsingDirectOrOpen(seed, index, adjustmentValue);
                }
            }
        } else {
            rank = extractRanking(seed);

            try {
                index = dirTable[codeTable[rank]];

                foundLocations = exactSearchingUsingDirectOrOpen(seed, index, adjustmentValue);
            } catch (exception& e) { }
        }
    }

    return foundLocations;
}

/**
* It gives the locations of the seed using approximate searching.
*
* @param seed - The seed that will be used for determining its locations in the reference genome.
* @param mode - It indicates what mode of hash-based indexer will be used.
* @param adjustmentValue - The value for the alignment of the locations found of the seed to the read.
* @return found locations of the seed
*/
vector<unsigned long long int> approximateSearchingPosition(string seed, string mode, int adjustmentValue) {
    vector<unsigned long long int> foundLocations;

    unsigned long long int rank;
    unsigned long long int index;

    if (mode.compare("min") == 0) {
        rank = getMinimizerRank(seed, q, w);
        vector<unsigned long long int> minRankLocations;
        minRankLocations = minimizers[rank];

        if (minRankLocations.size() > 0) {
            foundLocations = approximateSearchingUsingMinimizers(seed.substr(0, q), adjustmentValue, minRankLocations);
        }
    } else if (mode.compare("dir") == 0 || mode.compare("open") == 0) {
        if (mode.compare("dir") == 0) {
            rank = extractRanking(seed);

            if (rank >= 0) {
                index = dirTable[rank];

                if (index < posTable.size()) {
                    foundLocations = approximateSearchingUsingDirectOrOpen(seed, index, adjustmentValue);
                }
            }
        } else {
            rank = extractRanking(seed);

            try {
                index = dirTable[codeTable[rank]];

                foundLocations = approximateSearchingUsingDirectOrOpen(seed, index, adjustmentValue);
            } catch (exception& e) { }
        }
    }

    return foundLocations;
}