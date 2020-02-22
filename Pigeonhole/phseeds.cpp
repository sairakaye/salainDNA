//
// Created by saimanalili on 12/02/2020.
//

#include "phseeds.h"
#include "pigeonhole.h"
#include "minimizers.h"

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
