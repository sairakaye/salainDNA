//
// Created by saimanalili on 30/12/2019.
//

#include "pigeonhole.h"
#include "minimizers.h"

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

vector<AlignedReads *> processingRead(string read, int q, int m, int j, int k) {
    vector<AlignedReads *> reads;

    for (int i = 0; i < j; i++) {
        int startPosition = i * q;

        string seed = read.substr(startPosition, q);

        if (seed.length() < q) {
            //startPosition = (read.length() - 1) - q + 1;
            continue;
        }

        unsigned int rank = extractRanking(seed);

        vector<int> location = minimizers[rank];

        if (location.size() > 0) {
            for (int i2 = 0; i2 < location.size(); i2++) {
                if (seed.compare(refGenome.substr(location[i2], q)) == 0 && (location[i2] - startPosition) >= 0) {
                    AlignedReads *alignedRead = new AlignedReads;

                    alignedRead->startPos = location[i2] - startPosition;
                    alignedRead->read = read;
                    alignedRead->readFromGenome = refGenome.substr(location[i2] - startPosition, m);

                    reads.push_back(alignedRead);
                    //cout << read + " ";
                    //cout << refGenome.substr(location[i2] - startPosition, m) << endl;
                }
            }
        } else {
            continue;
        }
    }

    return reads;
}

//AlignedReads *processingRead(string read, int q, int m, int j) {
//    AlignedReads *reads = NULL;
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
//            for (int i2 = 0; i2 < location.size() && reads == NULL; i2++) {
//                if (read.substr(start, q).compare(refGenome.substr(location[i2], q)) == 0) {
//                    reads = new AlignedReads;
//
//                    reads->read = read;
//                    reads->readFromGenome = refGenome.substr(location[i2]-start, m);
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
//    return reads;
//}

void filterReads(string filename, int q, int m, int j, int k) {
    ofstream filterReads(filename + ".fpg");

    for (int i = 0; i < reads.size(); i++) {
        string forwardRead(reads[i]);

        vector<AlignedReads *> forward = processingRead(forwardRead, q, m, j, k);

        if (forward.size() > 0) {
            for (AlignedReads *alignedRead : forward) {
                filterReads << "[" + to_string(alignedRead->startPos) + "]" + alignedRead->read + "(Forward)" << endl;
                filterReads << to_string(alignedRead->startPos) + " " + alignedRead->read + " " + alignedRead->readFromGenome << endl;
                delete alignedRead;
            }
        } else {
            string reverseRead = reverseComplement(forwardRead);
            vector<AlignedReads *> reverse = processingRead(reverseRead, q, m, j, k);

            if (reverse.size() > 0) {
                for (AlignedReads *alignedRead : reverse) {
                    filterReads << "[" + to_string(alignedRead->startPos) + "]" + alignedRead->read + "(Reverse)" << endl;
                    filterReads << to_string(alignedRead->startPos) + " " + alignedRead->read + " " + alignedRead->readFromGenome << endl;
                    delete alignedRead;
                }
            } else {
                cout << "NOT INCLUDED! " + to_string(i) + " " + reads[i] << endl;
            }
        }

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
//                cout << "NOT INCLUDED! " + to_string(i) + " " + reads[i] << endl;
//            }
//        }
    }
}

//void parallelizeFilterReads(string filename, int q, int m, int j, int k) {
//    ofstream filterReads(filename + ".fpg");
//
//    #pragma omp parallel for
//    for (int i = 0; i < reads.size(); i++) {
//        string forwardRead(reads[i]);
//
//        AlignedReads *forward = processingRead(forwardRead, q, m, j, k);
//
//        if (forward != NULL) {
//            #pragma omp critical
//            filterReads << "[" + to_string(forward->startPos) + "]" + forward->read + "(Forward)" << endl;
//            filterReads << forward->read + " " + forward->readFromGenome << endl;
//            delete forward;
//        } else {
//            string reverseRead = reverseComplement(forwardRead);
//            AlignedReads *reverse = processingRead(reverseRead, q, m, j, k);
//
//            if (reverse != NULL) {
//                #pragma omp critical
//                filterReads << "[" + to_string(reverse->startPos) + "]" + reverse->read + "(Reverse)" << endl;
//                filterReads << reverse->read + " " + reverse->readFromGenome << endl;
//                delete reverse;
//            }
//        }
//    }
//}