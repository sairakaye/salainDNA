//
// Created by saimanalili on 12/02/2020.
//

#include "phreads.h"


/** DIFFERENT FUNCTIONS **/
vector<unsigned long long int> findingPositionUsingMinimizers(string windowSeed, string read, int windowLength, int q, int k, bool isForwardStrand) {
    vector<unsigned long long int> foundLocations;
    vector<unsigned long long int> location;
    unsigned int rank;

    if (windowSeed.length() == windowLength) {
        rank = getMinimizerRank(windowSeed, q, windowLength);
        location = minimizers[rank];

        if (location.size() > 0) {
            for (int i2 = 0; i2 < location.size(); i2++) {
                if (windowSeed.substr(0, q).compare(refGenome.substr(location[i2], q)) == 0) {
                    foundLocations.push_back(location[i2] - (windowLength * k));
                    readsMap[read].push_back(location[i2] - (q * k));

                    if (isForwardStrand) {
                        forwardFound.push_back(location[i2] - (q * k));
                    } else {
                        reverseFound.push_back(location[i2] - (q * k));
                    }
                }
            }
        }
    } else if (windowSeed.length() == q) {
        unsigned long long int rankHashValue  = extractRanking(windowSeed);
        location = minimizers[rankHashValue];

        if (location.size() > 0) {
            for (int i2 = 0; i2 < location.size(); i2++) {
                if (windowSeed.compare(refGenome.substr(location[i2], q)) == 0) {
                    foundLocations.push_back(location[i2] - (q * k));
                    readsMap[read].push_back(location[i2] - (q * k));

                    if (isForwardStrand) {
                        forwardFound.push_back(location[i2] - (q * k));
                    } else {
                        reverseFound.push_back(location[i2] - (q * k));
                    }
                }
            }
        }
    }

    return foundLocations;
}

vector<unsigned long long int> findingPositionUsingOpenAddr(string windowSeed, string read, int windowLength, int q, int k, bool isForwardStrand) {
    vector<unsigned long long int> foundLocations;
    vector<unsigned long long int> location;
    unsigned long long int rank;

    if (windowSeed.length() >= windowLength) {
        rank = extractRanking(windowSeed);

        try {
            //long long index2 = codeTable[rank];
            unsigned long long int index = dirTable[codeTable[rank]];

            while (windowSeed.compare(refGenome.substr(posTable[index], q)) == 0) {
                foundLocations.push_back(posTable[index] -  (q * k));
                readsMap[read].push_back(posTable[index] - (q * k));

                if (isForwardStrand) {
                    forwardFound.push_back(posTable[index] - (q * k));
                } else {
                    reverseFound.push_back(posTable[index] - (q * k));
                }

                index++;
            }
        } catch (exception& e) {

        }
    }

    return foundLocations;
}

vector<unsigned long long int> findingPositionUsingDirAddr(string windowSeed, string read, int windowLength, int q, int k, bool isForwardStrand) {
    vector<unsigned long long int> foundLocations;
    vector<unsigned long long int> location;
    unsigned long long int rank;

    if (windowSeed.length() >= windowLength) {
        rank = extractRanking(windowSeed);

        if (rank >= 0) {
            unsigned long long int index = dirTable[rank];

            if (index < posTable.size()) {
                while (windowSeed.compare(refGenome.substr(posTable[index], q)) == 0) {
                    foundLocations.push_back(posTable[index] - (q * k));
                    //foundLocations.push_back(posTable[index]);
                    readsMap[read].push_back(posTable[index] - (q * k));

                    if (isForwardStrand) {
                        forwardFound.push_back(posTable[index] - (q * k));
                    } else {
                        reverseFound.push_back(posTable[index] - (q * k));
                    }

                    index++;
                }
            }
        }
    }

    return foundLocations;
}
/*** END ***/

void partitioningReadsToSeeds(string filename, string mode, int windowLength, int q) {
    //ofstream filterSeeds(filename + "_" + to_string(q) + ".fpg");
    numberSeeds = 0;
    numberAcceptedSeeds = 0;
    numberOfFoundSeeds = 0;
    numberReads = 0;
    numberAcceptedReads = 0;

    bool isForward = false;

    int m = seeds[0].length();

    int j = m / q;

    //#pragma omp parallel for
    for (int i = 0; i < seeds.size(); i++) {
        string forwardRead(seeds[i]);
        numberReads++;

        bool isAccepted = false;

        vector<string> forward;
        vector<string> reverse;

        for (int k = 0; k < j; k++) {
            numberSeeds++;
            int startPosition = k * q;

            string seed;
            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, windowLength);
            } else {
                seed = forwardRead.substr(startPosition, q);
            }

            if (windowLength >= seed.size()) {
                vector<unsigned long long int> forward;
                isForward = true;

                if (mode.compare("min") == 0) {
                    forward = findingPositionUsingMinimizers(seed, seeds[i], windowLength, q, k, isForward);
                } else if (mode.compare("dir") == 0) {
                    forward = findingPositionUsingDirAddr(seed, seeds[i], windowLength, q, k, isForward);
                } else if (mode.compare("open") == 0) {
                    forward = findingPositionUsingOpenAddr(seed, seeds[i], windowLength, q, k, isForward);
                }

                if (forward.size() > 0) {
//                    for (unsigned long long int pos : forward) {
//                        //filterSeeds << fs << endl;
//                        EdlibAlignResult result = edlibAlign(refGenome.substr(pos, m).c_str(), m, seeds[i].c_str(), m, edlibDefaultAlignConfig());
//
//                        if (result.status == EDLIB_STATUS_OK) {
//                            if (result.editDistance == 0) {
//                                exactFound.push_back(pos);
//                            }
//                        }
//
//                        edlibFreeAlignResult(result);
//
//                        numberOfFoundSeeds++;
//                    }

                    //filterSeeds << forwardSeed << endl;
                    numberAcceptedSeeds++;

                    if (!isAccepted) {
                        isAccepted = true;
                    }
                } else {
                    cout << "REJECTED SEED! " + to_string(i + 1) + "-" + to_string(k) + " " + seed << endl;
                }
            } else {
                break;
            }
        }

        if (!isAccepted) {
            string reverseRead = reverseComplement(forwardRead);

            for (int k = 0; k < j; k++) {
                int startPosition = k * q;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, windowLength);
                } else {
                    seed = reverseRead.substr(startPosition, q);
                }

                if (windowLength >= seed.size()) {
                    vector<unsigned long long int> reverse;
                    isForward = false;

                    if (mode.compare("min") == 0) {
                        reverse = findingPositionUsingMinimizers(seed, seeds[i], windowLength, q, k, isForward);
                    } else if (mode.compare("dir") == 0) {
                        reverse = findingPositionUsingDirAddr(seed, seeds[i], windowLength, q, k, isForward);
                    } else if (mode.compare("open") == 0) {
                        reverse = findingPositionUsingOpenAddr(seed, seeds[i], windowLength, q, k, isForward);
                    }

                    if (reverse.size() > 0) {
//                        for (unsigned long long int pos : reverse) {
//                            EdlibAlignResult result = edlibAlign(refGenome.substr(pos, m).c_str(), m, seeds[i].c_str(), m, edlibDefaultAlignConfig());
//
//                            if (result.status == EDLIB_STATUS_OK) {
//                                if (result.editDistance == 0) {
//                                    exactFound.push_back(pos);
//                                }
//                            }
//
//                            edlibFreeAlignResult(result);
//
//                            numberOfFoundSeeds++;
//                        }

                        //filterSeeds << forwardSeed << endl;
                        numberAcceptedSeeds++;

                        if (!isAccepted) {
                            isAccepted = true;
                        }
                    } else {
                        cout << "REJECTED SEED! " + to_string(i + 1) + "-" + to_string(k) + " " + seed << endl;
                    }
                } else {
                    break;
                }
            }
        }

        if (isAccepted) {
            numberAcceptedReads++;
        }
    }

    cout << "Number of partitions per read: " + to_string(j) << endl;
    cout << "Total reads: " + to_string(numberReads) << endl;
    //cout << "Total accepted reads: " + to_string(numberAcceptedReads) << endl;
    //cout << "Total seeds evaluated: " + to_string(numberSeeds) << endl;
    //cout << "Total accepted seeds: evaluated " + to_string(numberAcceptedSeeds) << endl;
    //cout << "Total found seeds from the genome (non-unique): " + to_string(numberOfFoundSeeds) << endl;
}


/** VERSION TWO **/
void partitioningReadsToSeedsExit(string filename, string mode, int windowLength, int q) {
    //ofstream filterSeeds(filename + "_" + to_string(q) + ".fpg");
    numberSeeds = 0;
    numberAcceptedSeeds = 0;
    numberOfFoundSeeds = 0;
    numberReads = 0;
    numberAcceptedReads = 0;

    bool isForward = false;

    int m = seeds[0].length();

    int j = m / q;

    #pragma omp parallel for
    for (int i = 0; i < seeds.size(); i++) {
        string forwardRead(seeds[i]);
        numberReads++;

        bool isAccepted = false;

        vector<string> forward;
        vector<string> reverse;


        for (int k = 0; k < j; k++) {
            numberSeeds++;
            int startPosition = k * q;

            string seed;
            if (mode.compare("min") == 0) {
                seed = forwardRead.substr(startPosition, windowLength);
            } else {
                seed = forwardRead.substr(startPosition, q);
            }

            if (windowLength >= seed.size()) {
                vector<unsigned long long int> forward;

                if (mode.compare("min") == 0) {
                    forward = findingPositionUsingMinimizers(seed, seeds[i],  windowLength, q, k, isForward);
                } else if (mode.compare("dir") == 0) {
                    forward = findingPositionUsingDirAddr(seed, seeds[i], windowLength, q, k, isForward);
                } else if (mode.compare("open") == 0) {
                    forward = findingPositionUsingOpenAddr(seed, seeds[i], windowLength, q, k, isForward);
                }

                if (forward.size() > 0) {
//                    for (unsigned long long int pos : forward) {
//                        //filterSeeds << fs << endl;
//                        EdlibAlignResult result = edlibAlign(refGenome.substr(pos, m).c_str(), m, seeds[i].c_str(), m, edlibDefaultAlignConfig());
//
//                        if (result.status == EDLIB_STATUS_OK) {
//                            if (result.editDistance == 0) {
//
//                            }
//                        }
//
//                        edlibFreeAlignResult(result);
//                        numberOfFoundSeeds++;
//                    }

                    //filterSeeds << forwardSeed << endl;
                    numberAcceptedSeeds++;

                    if (!isAccepted) {
                        isAccepted = true;
                        break;
                    }
                } else {
                    cout << "REJECTED SEED! " + to_string(i + 1) + "-" + to_string(k) + " " + seed << endl;
                }
            } else {
                break;
            }
        }

        if (!isAccepted) {
            string reverseRead = reverseComplement(forwardRead);

            for (int k = 0; k < j; k++) {
                int startPosition = k * q;

                string seed;

                if (mode.compare("min") == 0) {
                    seed = reverseRead.substr(startPosition, windowLength);
                } else {
                    seed = reverseRead.substr(startPosition, q);
                }

                if (windowLength >= seed.size()) {
                    vector<unsigned long long int> reverse;


                    if (mode.compare("min") == 0) {
                        reverse = findingPositionUsingMinimizers(seed, seeds[i], windowLength, q, k, isForward);
                    } else if (mode.compare("dir") == 0) {
                        reverse = findingPositionUsingDirAddr(seed, seeds[i], windowLength, q, k, isForward);
                    } else if (mode.compare("open") == 0) {
                        reverse = findingPositionUsingOpenAddr(seed, seeds[i], windowLength, q, k, isForward);
                    }

                    if (reverse.size() > 0) {
//                        for (unsigned long long int pos : reverse) {
//                            EdlibAlignResult result = edlibAlign(refGenome.substr(pos, m).c_str(), m, seeds[i].c_str(), m, edlibDefaultAlignConfig());
//
//                            if (result.status == EDLIB_STATUS_OK) {
//                                if (result.editDistance == 0) {
//
//                                }
//                            }
//
//                            edlibFreeAlignResult(result);
//
//                            numberOfFoundSeeds++;
//                        }

                        //filterSeeds << forwardSeed << endl;
                        numberAcceptedSeeds++;

                        if (!isAccepted) {
                            isAccepted = true;
                            break;
                        }
                    } else {
                        cout << "REJECTED SEED! " + to_string(i + 1) + "-" + to_string(k) + " " + seed << endl;
                    }
                } else {
                    break;
                }
            }
        }

        if (isAccepted) {
            numberAcceptedReads++;
            //break;
        }
    }

    cout << "Number of partitions per read: " + to_string(j) << endl;
    cout << "Total reads: " + to_string(numberReads) << endl;
    cout << "Total accepted reads: " + to_string(numberAcceptedReads) << endl;
    cout << "Total seeds evaluated: " + to_string(numberSeeds) << endl;
    cout << "Total accepted seeds: evaluated " + to_string(numberAcceptedSeeds) << endl;
    cout << "Total found seeds from the genome (non-unique): " + to_string(numberOfFoundSeeds) << endl;
}
