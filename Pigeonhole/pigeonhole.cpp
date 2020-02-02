//
// Created by saimanalili on 30/12/2019.
//

#include <omp.h>
#include "pigeonhole.h"
#include "minimizers.h"

int numberSeeds = 0;
int numberAcceptedSeeds = 0;
unsigned int numberOfFoundSeeds = 0;
int numberReads = 0;
int numberAcceptedReads = 0;

string reverseComplement(string read) {
    string reverseRead = read;
    reverse(reverseRead.begin(), reverseRead.end());

    for (int i = 0; i < reverseRead.length(); i++) {
        switch (reverseRead[i]){
            case 'A':
                reverseRead[i] = 'T';
                break;
            case 'C':
                reverseRead[i] = 'G';
                break;
            case 'G':
                reverseRead[i] = 'C';
                break;
            case 'T':
                reverseRead[i] = 'A';
                break;
        }
    }

    return reverseRead;
}

//vector<AlignedReads *> processingRead(string read, int q, int m, int j, int k) {
//    vector<AlignedReads *> seeds;
//
//    for (int i = 0; i < j; i++) {
//        int startPosition = i * q;
//
//        string seed = read.substr(startPosition, q);
//
//        if (seed.length() < q) {
//            //startPosition = (read.length() - 1) - q + 1;
//            continue;
//        }
//
//        unsigned long rank = extractRanking(seed);
//
//        vector<unsigned long> location = minimizers[rank];
//
//        if (location.size() > 0) {
//            for (int i2 = 0; i2 < location.size(); i2++) {
//                if (seed.compare(refGenome.substr(location[i2], q)) == 0 && (location[i2] - startPosition) >= 0) {
//                    AlignedReads *alignedRead = new AlignedReads;
//
//                    alignedRead->startPos = location[i2] - startPosition;
//                    alignedRead->read = read;
//                    alignedRead->readFromGenome = refGenome.substr(location[i2] - startPosition, m);
//
//                    seeds.push_back(alignedRead);
//                    //cout << read + " ";
//                    //cout << refGenome.substr(location[i2] - startPosition, m) << endl;
//                }
//            }
//        } else {
//            continue;
//        }
//    }
//
//    return seeds;
//}

//vector<AlignedReads *> processingReadForOpen(string read, int q, int m, int j, int k) {
//    vector<AlignedReads *> seeds;
//
//    for (int i = 0; i < j; i++) {
//        int startPosition = i * q;
//
//        string seed = read.substr(startPosition, q);
//
//        if (seed.length() < q) {
//            //startPosition = (read.length() - 1) - q + 1;
//            continue;
//        }
//
//        /**
//         * For Open Addressing
//         */
//        int index = findIndex(read.substr(startPosition, q), indexing->codeTable);
//        int locationIndex = indexing->codeTable[index].rank;
//
//        while (refGenome.substr(indexing->posTable[locationIndex], q).compare(seed) == 0 &&
//               (indexing->posTable[locationIndex] - startPosition) >= 0) {
//            AlignedReads *alignedRead = new AlignedReads;
//            alignedRead->startPos = indexing->posTable[locationIndex];
//            alignedRead->read = read;
//            alignedRead->readFromGenome = refGenome.substr(indexing->posTable[locationIndex] - startPosition, m);
//
//            seeds.push_back(alignedRead);
//        }
//    }
//
//    return seeds;
//}

//vector<AlignedReads *> processingReadForDirect(string read, int q, int m, int j, int k) {
//    vector<AlignedReads *> seeds;
//
//    for (int i = 0; i < j; i++) {
//        int startPosition = i * q;
//
//        string seed = read.substr(startPosition, q);
//
//        if (seed.length() < q) {
//            //startPosition = (read.length() - 1) - q + 1;
//            continue;
//        }
//
//        /**
//         * For Direct Addressing
//         */
//        int rank = getRanking(read.substr(startPosition, q), q);
//        int locationIndex = directIndexing->dirTable[rank];
//
//        while (refGenome.substr(indexing->posTable[locationIndex], q).compare(seed) == 0 &&
//                (indexing->posTable[locationIndex] - startPosition) >= 0) {
//            AlignedReads *alignedRead = new AlignedReads;
//            alignedRead->startPos = indexing->posTable[locationIndex];
//            alignedRead->read = read;
//            alignedRead->readFromGenome = refGenome.substr(indexing->posTable[locationIndex] - startPosition, m);
//
//            seeds.push_back(alignedRead);
//        }
//    }
//
//    return seeds;
//}

//AlignedReads *processingRead(string read, int q, int m, int j) {
//    AlignedReads *seeds = NULL;
//
//    for (int i = 0; i < j; i++) {
//        int start = i * q;
//        int end = min(((i + 1) * q), (int) read.length());
//
//        if (end >= read.length()) {
//            end = read.length() - 1;
//            start = end - q + 1;
//        }
//
//        int rank = extractRanking(read.substr(start, q));
//
//        vector<int> location = minimizers[rank];
//
//        if (location.size() > 0) {
//            for (int i2 = 0; i2 < location.size() && seeds == NULL; i2++) {
//                if (read.substr(start, q).compare(refGenome.substr(location[i2], q)) == 0) {
//                    seeds = new AlignedReads;
//
//                    seeds->read = read;
//                    seeds->readFromGenome = refGenome.substr(location[i2]-start, m);
//
//                    cout << read + " ";
//                    cout << refGenome.substr(location[i2]-start, m) << endl;
//                }
//            }
//        } else {
//            continue;
//        }
//    }
//
//    return seeds;
//}

vector<string> processingAcceptingSeedUsingMinimizers(string windowSeed, int windowLength, int q, bool isForwardStrand) {
    vector<string> acceptedSeeds;
    vector<unsigned long int> location;
    unsigned int rank;

    if (windowSeed.length() >= windowLength) {
        rank = getMinimizerRank(windowSeed, q, windowLength);
        location = minimizers[rank];

        if (location.size() > 0) {
            for (int i2 = 0; i2 < location.size(); i2++) {
                if (windowSeed.compare(refGenome.substr(location[i2], windowLength)) == 0) {
                    acceptedSeeds.push_back(windowSeed);
                    //foundLocations.push_back(location[i2]);

                    if (isForwardStrand) {
                        forwardFound.push_back(location[i2]);
                    } else {
                        reverseFound.push_back(location[i2]);
                    }
                }
            }
        }
    } else if (windowSeed.length() == q) {
        unsigned long int rankHashValue  = extractRanking(windowSeed);
        location = minimizers[rankHashValue];

        if (location.size() > 0) {
            for (int i2 = 0; i2 < location.size(); i2++) {
                if (windowSeed.compare(refGenome.substr(location[i2], q)) == 0) {
                    acceptedSeeds.push_back(windowSeed);
                    //foundLocations.push_back(location[i2]);

                    if (isForwardStrand) {
                        forwardFound.push_back(location[i2]);
                    } else {
                        reverseFound.push_back(location[i2]);
                    }
                }
            }
        }
    }

    return acceptedSeeds;
}

vector<string> processingAcceptingSeedUsingOpenAddr(string windowSeed, int windowLength, int q, bool isForwardStrand) {
    vector<string> acceptedSeeds;
    vector<unsigned long int> location;
    unsigned long int rank;

    if (windowSeed.length() >= windowLength) {
        rank = extractRanking(windowSeed);

        try {
            //long long index2 = codeTable[rank];
            unsigned long int index = dirTable[codeTable[rank]];

            while (windowSeed.compare(refGenome.substr(posTable[index], q)) == 0) {
                acceptedSeeds.push_back(refGenome.substr(posTable[index], q));

                if (isForwardStrand) {
                    forwardFound.push_back(posTable[index]);
                } else {
                    reverseFound.push_back(posTable[index]);
                }

                index++;
            }
        } catch (exception& e) {

        }
    }

    return acceptedSeeds;
}

vector<string> processingAcceptingSeedUsingDirAddr(string windowSeed, int windowLength, int q, bool isForwardStrand) {
    vector<string> acceptedSeeds;
    vector<unsigned long int> location;
    unsigned long int rank;

    if (windowSeed.length() >= windowLength) {
        rank = extractRanking(windowSeed);

        if (rank >= 0) {
            unsigned long int index = dirTable[rank];

            if (index < posTable.size()) {
                while (windowSeed.compare(refGenome.substr(posTable[index], q)) == 0) {
                    acceptedSeeds.push_back(refGenome.substr(posTable[index], q));
                    //foundLocations.push_back(posTable[index]);

                    if (isForwardStrand) {
                        forwardFound.push_back(posTable[index]);
                    } else {
                        reverseFound.push_back(posTable[index]);
                    }

                    index++;
                }
            }
        }
    }

    return acceptedSeeds;
}

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
                vector<string> forward;
                mode = "forward";

                if (mode.compare("min") == 0) {
                    forward = processingAcceptingSeedUsingMinimizers(seed, windowLength, q, isForward);
                } else if (mode.compare("dir") == 0) {
                    forward = processingAcceptingSeedUsingDirAddr(seed, windowLength, q, isForward);
                } else if (mode.compare("open") == 0) {
                    forward = processingAcceptingSeedUsingOpenAddr(seed, windowLength, q, isForward);
                }

                if (forward.size() > 0) {
                    for (string fs : forward) {
                        //filterSeeds << fs << endl;
                        numberOfFoundSeeds++;
                    }

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
                    vector<string> reverse;


                    if (mode.compare("min") == 0) {
                        reverse = processingAcceptingSeedUsingMinimizers(seed, windowLength, q, isForward);
                    } else if (mode.compare("dir") == 0) {
                        reverse = processingAcceptingSeedUsingDirAddr(seed, windowLength, q, isForward);
                    } else if (mode.compare("open") == 0) {
                        reverse = processingAcceptingSeedUsingOpenAddr(seed, windowLength, q, isForward);
                    }

                    if (reverse.size() > 0) {
                        for (string fs : reverse) {
                            //filterSeeds << fs << endl;
                            numberOfFoundSeeds++;
                        }

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
    cout << "Total accepted reads: " + to_string(numberAcceptedReads) << endl;
    cout << "Total seeds evaluated: " + to_string(numberSeeds) << endl;
    cout << "Total accepted seeds: evaluated " + to_string(numberAcceptedSeeds) << endl;
    cout << "Total found seeds from the genome (non-unique): " + to_string(numberOfFoundSeeds) << endl;
}


void selectingSeeds(string filename, string mode, int windowLength, int q) {
    //ofstream filterSeeds(filename + "_" + to_string(q) + ".fpg");
    numberSeeds = 0;
    numberAcceptedSeeds = 0;
    numberOfFoundSeeds = 0;

    bool isForward = false;

    for (int i = 0; i < seeds.size(); i++) {
        string forwardSeed(seeds[i]);
        numberSeeds++;

        vector<string> forward;
        isForward = true;

        if (mode.compare("min") == 0) {
            forward = processingAcceptingSeedUsingMinimizers(forwardSeed, windowLength, q, isForward);
        } else if (mode.compare("dir") == 0) {
            forward = processingAcceptingSeedUsingDirAddr(forwardSeed, windowLength, q, isForward);
        } else if (mode.compare("open") == 0) {
            forward = processingAcceptingSeedUsingOpenAddr(forwardSeed, windowLength, q, isForward);
        }

        if (forward.size() > 0) {
            for (string fs : forward) {
                //filterSeeds << fs << endl;
                numberOfFoundSeeds++;
            }

            //filterSeeds << forwardSeed << endl;
            numberAcceptedSeeds++;
            forward.clear();
        } else {
            string reverseSeed = reverseComplement(forwardSeed);
            vector<string> reverse;
            isForward = false;

            if (mode.compare("min") == 0) {
                reverse = processingAcceptingSeedUsingMinimizers(reverseSeed, windowLength, q, isForward);
            } else if (mode.compare("dir") == 0) {
                reverse = processingAcceptingSeedUsingDirAddr(reverseSeed, windowLength, q, isForward);
            } else if (mode.compare("open") == 0) {
                reverse = processingAcceptingSeedUsingOpenAddr(reverseSeed, windowLength, q, isForward);
            }

            if (reverse.size() > 0) {
                for (string rs : reverse) {
                    //filterSeeds << rs << endl;
                    numberOfFoundSeeds++;
                }

                //filterSeeds << reverseSeed << endl;
                numberAcceptedSeeds++;
                reverse.clear();
            } else {
                cout << "REJECTED! " + to_string(i + 1) + " " + seeds[i] << endl;
            }

        }
    }

    cout << "Total seeds evaluated: " + to_string(numberSeeds) << endl;
    cout << "Total accepted seeds evaluated: " + to_string(numberAcceptedSeeds) << endl;
    cout << "Total found seeds from the genome (non-unique): " + to_string(numberOfFoundSeeds) << endl;
}



//void filterReads(string filename, int q, int m, int j, int k) {
//    ofstream filterReads(filename + "_" + to_string(q) +".fpg");
//
//    int rejectCount = 0;
//
//    for (int i = 0; i < seeds.size(); i++) {
//        string forwardRead(seeds[i]);
//
//        vector<AlignedReads *> forward = processingRead(forwardRead, q, m, j, k);
//        //vector<AlignedReads *> forward = processingReadForDirect(forwardRead, q, m, j, k);
//
//        if (forward.size() > 0) {
//            for (AlignedReads *alignedRead : forward) {
//                filterReads << "> Location: " + to_string(alignedRead->startPos) + " (Forward)" << endl;
//                filterReads << alignedRead->read + " " + alignedRead->readFromGenome << endl;
//                delete alignedRead;
//            }
//        } else {
//            string reverseRead = reverseComplement(forwardRead);
//            vector<AlignedReads *> reverse = processingRead(reverseRead, q, m, j, k);
//            //vector<AlignedReads *> reverse = processingReadForDirect(reverseRead, q, m, j, k);
//
//            if (reverse.size() > 0) {
//                for (AlignedReads *alignedRead : reverse) {
//                    filterReads << "Location: " + to_string(alignedRead->startPos) + " (Reverse)" << endl;
//                    filterReads << alignedRead->read + " " + alignedRead->readFromGenome << endl;
//                    delete alignedRead;
//                }
//            } else {
//                rejectCount++;
//                cout << "NOT INCLUDED! " + to_string(i) + " " + seeds[i] << endl;
//            }
//        }

//        AlignedReads *forward = processingRead(forwardRead, q, m, j, k);
//
//        if (forward != NULL) {
//            //filterReads << "[" + to_string(forward->startPos) + "]" + forward->read + "(Forward)" << endl;
//            filterReads << to_string(forward->startPos) + " " + forward->read + " " + forward->readFromGenome << endl;
//            delete forward;
//        } else {
//            string reverseRead = reverseComplement(forwardRead);
//            AlignedReads *reverse = processingRead(reverseRead, q, m, j, k);
//
//            if (reverse != NULL) {
//                //filterReads << "[" + to_string(reverse->startPos) + "]" + reverse->read + "(Reverse)" << endl;
//                filterReads << to_string(reverse->startPos) + " " + reverse->read + " " + reverse->readFromGenome << endl;
//                delete reverse;
//            } else {
//                cout << "NOT INCLUDED! " + to_string(i) + " " + seeds[i] << endl;
//            }
//        }
//    }
//}

//void parallelizeFilterReads(string filename, int q, int m, int j, int k) {
//    ofstream filterReads(filename + ".fpg");
//
//    omp_set_dynamic(0);
//    #pragma omp parallel for num_threads(4)
//    for (int i = 0; i < seeds.size(); i++) {
//        string forwardRead(seeds[i]);
//
//        vector<AlignedReads *> forward = processingRead(forwardRead, q, m, j, k);
//
//        if (forward.size() > 0) {
//            #pragma omp critical
//            for (AlignedReads *alignedRead : forward) {
//                filterReads << "[" + to_string(alignedRead->startPos) + "]" + alignedRead->read + "(Forward)" << endl;
//                filterReads << to_string(alignedRead->startPos) + " " + alignedRead->read + " " + alignedRead->readFromGenome << endl;
//                delete alignedRead;
//            }
//        } else {
//            string reverseRead = reverseComplement(forwardRead);
//            vector<AlignedReads *> reverse = processingRead(reverseRead, q, m, j, k);
//
//            if (reverse.size() > 0) {
//                #pragma omp critical
//                for (AlignedReads *alignedRead : reverse) {
//                    filterReads << "[" + to_string(alignedRead->startPos) + "]" + alignedRead->read + "(Reverse)" << endl;
//                    filterReads << to_string(alignedRead->startPos) + " " + alignedRead->read + " " + alignedRead->readFromGenome << endl;
//                    delete alignedRead;
//                }
//            } else {
//                cout << "NOT INCLUDED! " + to_string(i) + " " + seeds[i] << endl;
//            }
//        }

//        AlignedReads *forward = processingRead(forwardRead, q, m, j, k);
//
//        if (forward != NULL) {
//            //filterReads << "[" + to_string(forward->startPos) + "]" + forward->read + "(Forward)" << endl;
//            filterReads << to_string(forward->startPos) + " " + forward->read + " " + forward->readFromGenome << endl;
//            delete forward;
//        } else {
//            string reverseRead = reverseComplement(forwardRead);
//            AlignedReads *reverse = processingRead(reverseRead, q, m, j, k);
//
//            if (reverse != NULL) {
//                //filterReads << "[" + to_string(reverse->startPos) + "]" + reverse->read + "(Reverse)" << endl;
//                filterReads << to_string(reverse->startPos) + " " + reverse->read + " " + reverse->readFromGenome << endl;
//                delete reverse;
//            } else {
//                cout << "NOT INCLUDED! " + to_string(i) + " " + seeds[i] << endl;
//            }
//        }
//    }
//}