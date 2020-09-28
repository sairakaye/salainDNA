/**
 * @file bitmatrix.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the Bit Matrix Filter Process implementation.
 */

#include "bitmatrix.h"

int truePos; // The number  of True Positive locations filtered.
int trueNeg; // The number of True Negative locations filtered.

unsigned int alignmentNeeded = 0; // The number of accepted filtered locations.
unsigned int notNeeded = 0; // The number of rejected filtered locations.

/**
 * It counts the number of ones.
 *
 * @param toCount - A list of numbers to count.
 * @param E - The value of error threshold.
 * @return count of ones
 */
int countOnes(vector<int> toCount, unsigned int E) {
    int ctr = 0;

    if (toCount.empty()) {
        return E + 1;
    } else {
        for (int i = 0; i < toCount.size(); i++) {
            if (toCount[i] == 1) {
                ctr++;
            }
        }
    }
    return ctr;
}

/**
 * It does the bit matrix algorithm.
 *
 * @param E - The value of error threshold.
 * @param read - The read to be used for checking.
 * @param genome - The data of the reference genome.
 * @return the processed read with values 0 or 1.
 */
vector<int> bitMatrixAlgorithm(int E, string read, string genome) {

    vector<int> locations(2 * E + 1);
    int windowSize = 4;
    vector<int> finalVector(m);
    int m = genome.size();
    int count = 0;
    int diagLoc = 0;

    vector<int> NMap((2 * E + 1) * m + (2 * E + 1));
    vector<int> bestDiagonal = {-1, -1, -1, -1};
    vector<int> movingLocation;

    //NEIGHBORHOOD MAP
    int p = 0;
    for (int i = 0; i < 2 * E + 1; i++) {
        locations[i] = p;
        NMap[p] = 2;
        p += m + 1;
    }
    movingLocation = locations;

    int jE;
    for (int i = 0; i < m; i++) {
        jE = 0;
        for (int j = i - E; j < i + E + 1; j++) {
            if (jE <= 2 * E) {
                if ((genome[j] == 'A' || genome[j] == 'C' || genome[j] == 'G' || genome[j] == 'T') &&
                    read[i] == genome[j]) {
                    NMap[movingLocation[jE] + 1] = 0;

                    if (jE == E + 1) {
                        finalVector[diagLoc] = 0;
                        diagLoc++;
                    }
                    movingLocation[jE] += 1;
                    jE++;
                } else if ((genome[j] == 'A' || genome[j] == 'C' || genome[j] == 'G' || genome[j] == 'T') &&
                           read[i] != genome[j]) {
                    NMap[movingLocation[jE] + 1] = 1;

                    if (jE == E + 1) {
                        finalVector[diagLoc] = 1;
                        diagLoc++;
                        count++;
                    }

                    movingLocation[jE] += 1;
                    jE++;
                } else {
                    movingLocation[jE] += 1;
                    jE++;
                }
            }
        }
    }

    if (E == 0 || count <= E) {
        finalVector = NMap;
    } else {
        // END OF NEIGHBORHOOD MAP
        // SLIDING WINDOW
        int ctr = 0;
        for (int w = 0; w < m; w++) {
            if (w + windowSize > m) {
                windowSize = m - w;
            }

            // BEST DIAGONAL
            int bestZeroCount = 0;
            int currZeroCount = 0;
            int index = 0;

            for (int i = 0; i < E * 2 + 1; i++) {

                for (int j = 0; j < windowSize; j++) {
                    if (NMap[locations[i] + 1 + j + ctr] == 0) {
                        currZeroCount++;
                    }
                }
                if (currZeroCount >= bestZeroCount) {
                    bestZeroCount = currZeroCount;
                    index = i;
                }
                currZeroCount = 0;
            }
            //END OF BEST DIAGONAL

            ctr++;
            int l = 0;
            for (int i = w; i < windowSize + w; i++) {
                if (i < m) {
                    finalVector[i] = NMap[locations[index] + 1 + l + ctr];
                    l++;

                }
            }
        }
    }
    //END OF SLIDING WINDOW

    return finalVector;
}

/**
 * It starts the bit matrix filter process.
 *
 */
void bitMatrixFilterProcess() {
    alignmentNeeded = 0;
    notNeeded = 0;

    auto start = std::chrono::high_resolution_clock::now();

    int i;
    #pragma omp parallel for reduction(+:notNeeded, alignmentNeeded)
    for (i = 0; i < reads.size(); i++) {
        string read = reads[i].readData;

        if (reads[i].forwardLocations.size() > 0) {
            int j;
            for (j = 0; j < reads[i].forwardLocations.size(); j++) {
                if (refGenome.genomeData.substr(reads[i].forwardLocations[j], m).size() == m) {
                    if (countOnes(
                            bitMatrixAlgorithm(e, read, refGenome.genomeData.substr(reads[i].forwardLocations[j], m)), e) <= e) {
                        alignmentNeeded++;
                    } else {
                        notNeeded++;
                        #pragma omp critical
                        reads[i].forwardLocations.erase(reads[i].forwardLocations.begin() + j);
                        j--;
                    }
                } else {
                    notNeeded++;
                    #pragma omp critical
                    reads[i].forwardLocations.erase(reads[i].forwardLocations.begin() + j);
                    j--;
                }
            }
        }

        if (reads[i].reverseLocations.size() > 0) {
            string reverseRead = reverseComplement(read);

            int j;
            for (j = 0; j < reads[i].reverseLocations.size(); j++) {
                if (refGenome.genomeData.substr(reads[i].reverseLocations[j], m).size() == m) {
                    if (countOnes(bitMatrixAlgorithm(e, reverseRead,
                        refGenome.genomeData.substr(reads[i].reverseLocations[j], m)), e) <= e) {
                        #pragma omp critical
                        alignmentNeeded++;
                    } else {
                        notNeeded++;
                        #pragma omp critical
                        reads[i].reverseLocations.erase(reads[i].reverseLocations.begin() + j);
                        j--;
                    }
                } else {
                    notNeeded++;
                    #pragma omp critical
                    reads[i].reverseLocations.erase(reads[i].reverseLocations.begin() + j);
                    j--;
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();

    chrono::duration<double> diff = end - start;
    bmRunTime = diff.count();

    numFilteredReadLocations = alignmentNeeded;

    cout << "Time taken by the Bit Matrix is: " << to_string(bmRunTime) << " sec" << endl;
    cout << "Number of needed locations: " << to_string(alignmentNeeded) << endl;
    cout << "Number of not needed locations: " << to_string(notNeeded) << endl << endl;
}