//
// Created by saimanalili on 30/12/2019.
//

#include "pigeonhole.h"
#include "indexing.h"
#include "minimizers.h"

using namespace std;

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

bool isReadAccepted(string read, RefGenome *refGenome, HashIndexing *hashIndexing, int e, int segmentLength) {
    for (int i = 0; i < e + 1; i++) {
        int start = i * segmentLength;
        int end = min(((i + 1) * segmentLength), (int) read.length());

        // If the partition is not given the even length.
        if (end >= read.length()) {
            end = read.length() - 1;
            start = end - segmentLength + 1;
        }

        /* THIS IS FOR THE INDEX */
        long rank = getRanking(read.substr(start, end), segmentLength);
        int rankPos = findIndex(hashIndexing->codeTable, hashIndexing->codeTableLength, rank);

        if (rankPos < 0)
            continue;

        int posTableIndex = hashIndexing->dirTable[rankPos];
        int startElement = hashIndexing->posTable[posTableIndex];

        char *temp = (char *)read.c_str();

        while (strcmp(substring(temp, start, segmentLength), substring(refGenome->genome, startElement, segmentLength)) == 0) {
            return true;

//            posTableIndex++;
//
//            if (posTableIndex < hashIndexing->posTableLength) {
//                startElement = hashIndexing->posTable[posTableIndex];
//            } else {
//                cout << substring(temp, start, segmentLength) << endl;
//                break;
//            }
        }

        /* END HERE */
    }

    return false;
}

// e here is the j hehe
// segmentLength = m
AlignedReads *isReadAcceptedWithMinimizers(string read, RefGenome *refGenome, map<int, vector<int>> minimizers, int e, int segmentLength) {
    AlignedReads *reads = NULL;

    for (int i = 0; i < e; i++) {
        int start = i * segmentLength;
        int end = min(((i + 1) * segmentLength), (int) read.length());

        // If the partition is not given the even length.
        if (end >= read.length()) {
            end = read.length() - 1;
            start = end - segmentLength + 1;
        }

        /* USE MINIMIZERRRRR!!! */
        //int rank = generateMinimizerHash(read.substr(start, end), 8, 8);
        int rank = extractRanking(read.substr(start, segmentLength));

        vector<int> location = minimizers[rank];

        if (location.size() > 0) {
            for (int j = 0; j < location.size(); j++) {
                char *temp = (char *)read.substr(start, segmentLength).c_str();

                if (strcmp(temp, substring(refGenome->genome, location[j], segmentLength)) == 0) {
                    reads = (AlignedReads *)malloc(sizeof(AlignedReads));

                    reads->read = (char *)read.c_str();
                    reads->readFromGenome = substring(refGenome->genome, location[j]-start, 120);
                    //locationToGenome.push_back(location[j]);
                    cout << "[" + read + "] ";
                    cout << "[" + string(substring(refGenome->genome, location[j]-start, 120)) + "] " << endl;
                }
            }
        } else {
            continue;
        }
    }

    //return locationToGenome;
    return reads;
}

void filterReadsWithMinimizers(string filename, ReadList *readList, map<int, vector<int>> minimizers, RefGenome *refGenome, int e, int segmentLength) {
    ofstream filterReads(filename + ".fpg");

    for (int i = 0; i < readList->size; i++) {
        string forwardRead(readList->reads[i]);

        //vector<int> toAcceptForward = isReadAcceptedWithMinimizers(forwardRead, refGenome, minimizers, e, segmentLength);
        //bool toAcceptForward = isReadAcceptedWithMinimizers(forwardRead, refGenome, minimizers, e, segmentLength);

        AlignedReads *toAcceptForward = toAcceptForward = isReadAcceptedWithMinimizers(forwardRead, refGenome, minimizers, e, segmentLength);

//        if (toAcceptForward.size() > 0) {
//            stringstream result;
//            copy(toAcceptForward.begin(), toAcceptForward.end(), ostream_iterator<int>(result, " "));
//
//            filterReads << forwardRead + " (Forward) [" + result.str() + "]" << endl;
//        } else {
//            string reverseRead = reverseComplement(forwardRead);
//            //vector<int> toAcceptBackward = isReadAcceptedWithMinimizers(reverseRead, refGenome, minimizers, e, segmentLength);
//            bool toAcceptBackward = isReadAcceptedWithMinimizers(reverseRead, refGenome, minimizers, e, segmentLength);
//
//            if (toAcceptBackward.size() > 0) {
//                stringstream result;
//                copy(toAcceptBackward.begin(), toAcceptBackward.end(), ostream_iterator<int>(result, " "));
//
//                filterReads << reverseRead + " (Reverse) [ " + result.str() + "]" << endl;
//            }
//        }

        if (toAcceptForward) {
            //filterReads << forwardRead + " (Forward)" << endl;
            filterReads << string(toAcceptForward->read) + " " + string(toAcceptForward->readFromGenome) << endl;
            free(toAcceptForward);
        } else {
            string reverseRead = reverseComplement(forwardRead);
            AlignedReads *toAcceptBackward = isReadAcceptedWithMinimizers(reverseRead, refGenome, minimizers, e, segmentLength);

            if (toAcceptBackward) {
                //filterReads << reverseRead + " (Reverse)" << endl;
                filterReads << string(toAcceptBackward->read) + " " + string(toAcceptBackward->readFromGenome) << endl;
                free(toAcceptBackward);
            }
        }
    }
}

void filterReads(string filename, ReadList *readList, HashIndexing *hashIndexing, RefGenome *refGenome, int e, int segmentLength)  {
    ofstream filterReads(filename + ".fpg");

    for (int i = 0; i < readList->size; i++) {
        string forwardRead(readList->reads[i]);
        string reverseRead = reverseComplement(forwardRead);

        bool toAcceptForward = isReadAccepted(forwardRead, refGenome, hashIndexing, e, segmentLength);
        bool toAcceptBackward = isReadAccepted(reverseRead, refGenome, hashIndexing, e, segmentLength);

        if (toAcceptForward) {
            filterReads << forwardRead + " (Forward)" << endl;
        } else if (toAcceptBackward) {
            filterReads << reverseRead + " (Reverse)" << endl;
        }
    }
}

void parallelizeFilterReads(string filename, ReadList *readList, HashIndexing *hashIndexing, RefGenome *refGenome, int e, int segmentLength) {
    ofstream filterReads(filename + ".fpg");

    #pragma omp parallel for
    for (int i = 0; i < readList->size; i++) {
        string compareRead(readList->reads[i]);

        bool toAcceptRead = isReadAccepted(compareRead, refGenome, hashIndexing, e, segmentLength);

        #pragma omp critical
        if (toAcceptRead) {
            filterReads << compareRead + " " + " is accepted!" << endl;
        } else {
            filterReads << compareRead + " " + " is rejected!" << endl;
        }
    }
}