#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>

#include "file_manager.h"
#include "indexing.h"
#include "pigeonhole.h"

using namespace std;

RefGenome *refGenome;
map<int, vector<int>> minimizers;
ReadList *readList;

int main(int argc, char *argv[]) {
    cout << "Processing..." << endl;

    string filepath = "sim_255.fa";

    //string filepath = "testing.fa";

    filepath.erase(remove(filepath.begin(), filepath.end(), '\"'), filepath.end());
    filepath.erase(remove(filepath.begin(), filepath.end(), '\''), filepath.end());

    refGenome = readGenome(filepath);
    int m = 120;
    int q = 8;
    int k = 2;

    string fileMinimizers = "min_255.txt";
    minimizers = getMinimizersFromFile(fileMinimizers);

    cout << "Getting minimizers done!" << endl;

    /*
    float loadFactor = 0.8;

    hashIndexing = new HashIndexing(refGenome->length, k, loadFactor);

    processIndexing(hashIndexing, refGenome);
    cout << "Indexing complete!" << endl;


    cout << "Contents of the Indexing!" << endl;

    cout << "Code Table" << endl;
    for (int i = 0; i < hashIndexing->codeTableLength; i++) {
        cout << "Index: ";
        cout << i;
        cout << " Rank: ";
        cout << hashIndexing->codeTable[i].rank << endl;
    }

    cout << "Directory Table" << endl;
    for (int i = 0; i < hashIndexing->dirTableLength; i++) {
        cout << hashIndexing->dirTable[i];
        cout << " ";
    }

    cout << "Position Table" << endl;

    for (int i = 0; i < hashIndexing->posTableLength; i++) {
        cout << "Index: ";
        cout << i;
        cout << " Position: ";
        cout << hashIndexing->posTable[i] << endl;
    }
    */

    /** TO DO HERE : PIGEONHOLE PROCESS */
    string readsFilename = "sim_255_120bp_10r_perf.fa";
    //string readsFilename = "testing_reads.fa";
    readList = readReads(readsFilename);

    int j = m / q;
    long start = clock();
    filterReadsWithMinimizers(readsFilename, readList, minimizers, refGenome, j, q);
    //parallelizeFilterReads(readsFilename, readList, hashIndexing, refGenome, e, k);
    long end = clock();
    cout << "Pigeonhole process done!" << endl;

    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by the pigeonhole process is : " << fixed
         << time_taken << " sec" << endl;

    /** CLEANING UP OF EVERYTHING **/
    if (refGenome != NULL) {
        free(refGenome->genome);
        free(refGenome);
    }

//    if (readList != NULL) {
//        for (int i = 0; i < readList->size; i++) {
//            free(readList->reads[i]);
//        }
//
//        free(readList);
//    }

    return 0;
}