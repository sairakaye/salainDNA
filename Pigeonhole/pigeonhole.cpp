//
// Created by saimanalili on 30/12/2019.
//

#include <omp.h>
#include "pigeonhole.h"
#include "edlib.h"

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

void verification(vector<unsigned long long int>& forwardFound, vector<unsigned long long int>& reverseFound, vector<unsigned long long int>& exactFound) {
    cout << "Removing duplicate locations of the read/seeds mapped..." << endl;
    sort(forwardFound.begin(), forwardFound.end());
    forwardFound.erase(unique(forwardFound.begin(), forwardFound.end()), forwardFound.end());

    sort(reverseFound.begin(), reverseFound.end());
    reverseFound.erase(unique(reverseFound.begin(), reverseFound.end()), reverseFound.end());

    cout << "Number of seeds found from forward (unique locations): " + to_string(forwardFound.size()) << endl;
    cout << "Number of seeds found from backward (unique locations): " + to_string(reverseFound.size()) << endl;

    vector<unsigned long long int> combined(forwardFound);
    combined.insert(combined.end(), reverseFound.begin(), reverseFound.end());
    sort(combined.begin(), combined.end());
    combined.erase(unique(combined.begin(), combined.end()), combined.end());

    cout << "Number of seed locations from forward and backward (including intersection): " + to_string(forwardFound.size() + reverseFound.size()) << endl;
    cout << "Number of seed locations accepted (from both): " + to_string(combined.size()) << endl;


    sort(exactFound.begin(), exactFound.end());
    exactFound.erase(unique(exactFound.begin(), exactFound.end()), exactFound.end());

    cout << "Number of locations found the exact reads: " + to_string(exactFound.size()) << endl;

    for (unsigned long long int pos : exactFound) {
        cout << to_string(pos) << endl;
    }
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