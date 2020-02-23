#include "common.h"
#include "phreads.h"
#include "phseeds.h"
#include "edlib.h"

string refGenome;
vector<string> seeds;
//vector<unsigned long long int> foundLocations;
map<unsigned long long int, vector<unsigned long long int>> minimizers;
map<long long, unsigned long long int> codeTable;
vector<unsigned long long int> dirTable;
vector<unsigned long long int> posTable;
vector<unsigned long long int> forwardFound;
vector<unsigned long long int> reverseFound;
vector<unsigned long long int> exactFound;
map<string, vector<unsigned long long int>> readsMap;

void initializeMinimizersFromFile(string filename, map<unsigned long long int, vector<unsigned long long int>>& minimizers) {
    //cout << "Processing minimizers..." << endl;
    minimizers = getMinimizersFromFile(filename);
    //cout << "Done!" << endl;
}

void initializeDirectAddressingFromFile(string filename, vector<unsigned long long int>& dirTable,
        vector<unsigned long long int>& posTable) {
    //cout << "Processing direct addressing..." << endl;
    getDirectAddressing(filename, dirTable, posTable);
    //cout << "Done!" << endl;
}

void initializingOpenAddressingFromFile(string filename,  map<long long, unsigned long long int>& codeTable,
        vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable) {
    //cout << "Processing open addressing..." << endl;
    getOpenAddressing(filename, codeTable, dirTable, posTable);
    //cout << "Done!" << endl;
}

int main(int argc, char *argv[]) {
    //int w = 12;
    //int windowLength = w + q - 1;
    //int m = 100;
    //int k = 2;
    //int j = m / q;
    //float loadFactor = 0.8;

    /****** EDIT THIS SECTION ONLY (REPORT IF THERE ARE BUGS) *******/

    /**
     * m - length of the read.
     * q - length of the seed (k-mer/q-gram)
     * k - errors
     * j - number of partitions in the read.
     */
    int q = 8;
    int windowLength = q + q - 1;

    /**
     * Genome (1MB, 2MB, 4MB, 8MB)
     * - Take note that this is CASE SENSITIVE
     */
     string genomeType = "chr04_NoN";

    /**
     * isRead - make it true to enable the reading of the reads with 100bp in length;
     */
    bool isRead = true;

    /**
     * Indexing modes:
     * min - minimizer
     * dir - direct addressing
     * open - open addressing
     */
    string mode = "min";

    /**
     * Type of reads (only input these)
     * - perfect (can be found in the reference)
     * - random (cannot be found in the reference)
     * - mixed (50% perfect, 50% random)
     *
     */
    string readType = "perfect";

    /**
     * Number of reads or seeds. Only 1 and 100000.
     */
     int numR = 100;

    /*******               END OF SECTION                *******/

    //    filepath.erase(remove(filepath.begin(), filepath.end(), '\"'), filepath.end());
    //    filepath.erase(remove(filepath.begin(), filepath.end(), '\''), filepath.end());

    //string filepath = "./genome/ndna_" + genomeType + ".fa";
    string filepath = "./genome/" + genomeType + ".fa";
    cout << "Reading the reference genome... " << endl << filepath << endl;
    refGenome = readGenome(filepath);
    //refGenome = "TTATCTCTTA";
    cout << "Done!" << endl << endl;

    string readsFilename;

    if (isRead) {
        //readsFilename = "./reads/" + readType + "/" + readType + "_" + to_string(numR) + "r_" + "100bp.fa";
        readsFilename = "/home/saimanalili/multicore-rm/Pigeonhole/Reads - real data set/chr04_NoN_1_150R.fa";
    } else {
        if (mode.compare("dir") == 0 || mode.compare("open") == 0) {
            //readsFilename = "./seeds/for_direct_n_open/" + readType + "/" + readType + "_" + to_string(numR) + "r_" + to_string(q) + "bp.fa";
            readsFilename = "./seeds/for_direct_n_open/" + readType + "/" + readType + "_" + "chr04_" + to_string(numR) + "r_" + to_string(q) + "bp.fa";
        } else if (mode.compare("min") == 0) {
            //readsFilename = "./seeds/for_minimizers/" + readType + "/" + readType + "_" + to_string(numR) + "r_" + to_string(windowLength) + "bp.fa";
            readsFilename = "./seeds/for_minimizers/" + readType + "/" + readType + "_" + "chr04_" + to_string(numR) + "r_" + to_string(windowLength) + "bp.fa";
        } else {
            cout << "INVALID!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    cout << "Reading the reads/seeds... " << endl << readsFilename << endl;
    seeds = readReads(readsFilename);
    cout << "Done!" << endl << endl;

//    string fileMinimizers = "./min/min_" + genomeType + "_" + to_string(q) + ".txt";
//    string fileOpenAddr = "./open/open_" + genomeType + "_" + to_string(q) + ".txt";
//    string fileDirAddr = "./dir/dir_" + genomeType + "_" + to_string(q) + ".txt";
    string fileMinimizers = "./min/min_chr04_" + to_string(q) + ".txt";
    string fileOpenAddr = "./open/open_chr04_" + to_string(q) + ".txt";
    string fileDirAddr = "./dir/dir_chr04_" + to_string(q) + ".txt";

    if (mode.compare("min") == 0) {
        cout << "Initializing minimizers... " << endl << fileMinimizers << endl;
        initializeMinimizersFromFile(fileMinimizers, minimizers);
        cout << "Done!" << endl << endl;

        cout << "Doing pigeonhole process... " << endl;
        long start = clock();
        if (isRead) {
            partitioningReadsToSeeds(readsFilename, mode, windowLength, q);
        } else {
            selectingSeeds(readsFilename, mode, windowLength, q);
        }
        long end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        cout << "Time taken by the pigeonhole process is : " << fixed
             << time_taken << " sec" << endl;
        cout << "Done!" << endl << endl;
    } else if (mode.compare("dir") == 0) {
        cout << "Initializing direct addressing... " << endl << fileDirAddr << endl;
        initializeDirectAddressingFromFile(fileDirAddr, dirTable, posTable);
        cout << "Done!" << endl << endl;

        cout << "Doing pigeonhole process..." << endl;
        long start = clock();
        if (isRead) {
            partitioningReadsToSeeds(readsFilename, mode, q, q);
            //partitioningReadsToSeedsExit(readsFilename, mode, q, q);
        } else {
            selectingSeeds(readsFilename, mode, q, q);
        }
        long end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        cout << "Time taken by the pigeonhole process is : " << fixed
             << time_taken << " sec" << endl;
        cout << "Done!" << endl << endl;
    } else if (mode.compare("open") == 0) {
        cout << "Initializing open addressing... " << endl << fileOpenAddr << endl;
        initializingOpenAddressingFromFile(fileOpenAddr, codeTable, dirTable, posTable);
        cout << "Done!" << endl << endl;

        cout << "Doing pigeonhole process..." << endl;
        long start = clock();
        if (isRead) {
            partitioningReadsToSeeds(readsFilename, mode, q, q);
        } else {
            selectingSeeds(readsFilename, mode, q, q);
        }
        long end = clock();

        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        cout << "Time taken by the pigeonhole process is : " << fixed
             << time_taken << " sec" << endl;
        cout << "Done!" << endl << endl;

    } else {
        cout << "INVALID!" << endl;
        exit(EXIT_FAILURE);
    }

    for (pair<string, vector<unsigned long long int>> readPair : readsMap) {
        for (unsigned long long int location : readPair.second) {
            EdlibAlignResult result = edlibAlign(refGenome.substr(location, readPair.first.length()).c_str(), readPair.first.length(), readPair.first.c_str(), readPair.first.length(), edlibDefaultAlignConfig());

            if (result.status == EDLIB_STATUS_OK) {
                if (result.editDistance == 0) {
                    exactFound.push_back(location);
                }

                edlibFreeAlignResult(result);
            }
        }
    }

    verification(forwardFound, reverseFound, exactFound);

    return 0;
}