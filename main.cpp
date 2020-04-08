#include <omp.h>
#include "common.h"
#include "pigeonhole.h"

#include "bitmatrix.h"
#include <chrono>

using namespace std::chrono;

string refGenome;
vector<string> reads;
map<string, vector<unsigned long long int>> forwardReadsMap;
map<string, vector<unsigned long long int>> reverseReadsMap;

map<unsigned long long int, vector<unsigned long long int>> minimizers;
map<long long, unsigned long long int> codeTable;
vector<unsigned long long int> dirTable;
vector<unsigned long long int> posTable;

string temp_comp;
ofstream outputPossibleReadsFile;
ofstream outputLocationsFile;
ofstream infoFile;

Counters *counter;

string mode;
unsigned int q;
unsigned int w;
unsigned int m;
unsigned int e;

int main(int argc, char *argv[]) {
    string genomeFileName;
    string readsFileName;
    string indexFileName;

    q = 8;
    e = 0;
    mode = "min";
    string searchMode = "all";
    /**
     * TO DO: Argument checking... (if missing arguments etc.)
     * ALSO TO DO: accept paths with spaces.
     */
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]) == "-q") {
            q = atoi(argv[i + 1]);
        } else if (string(argv[i]) == "-g") {
            genomeFileName = string(argv[i + 1]);
        } else if (string(argv[i]) == "-ir") {
            readsFileName = string(argv[i + 1]);
        } else if (string(argv[i]) == "-i") {
            indexFileName = string(argv[i + 1]);
        } else if (string(argv[i]) == "-m") {
            mode = string(argv[i + 1]);
        } else if (string(argv[i]) == "-temp") {
            temp_comp = string(argv[i + 1]);
        } else if (string(argv[i]) == "-e") {
            e = atoi(argv[i + 1]);
        } else if (string(argv[i]) == "-s") {
            searchMode = string(argv[i + 1]);
        }
    }

    cout << "Reading the reference genome... " << endl << genomeFileName << endl << endl;
    refGenome = readGenomeFile(genomeFileName);
    cout << "Reading the reads... " << endl << readsFileName << endl << endl;
    reads = readReadsFile(readsFileName);

    w = q + q - 1;
    m = reads[0].size();

    cout << "Reading the indexing... " << endl << indexFileName << endl << endl;
    if (mode.compare("min") == 0) {
        minimizers = getMinimizers(indexFileName);
    } else if (mode.compare("dir") == 0) {
        getDirectAddressing(indexFileName, dirTable, posTable);
    } else if (mode.compare("open") == 0) {
        getOpenAddressing(indexFileName, codeTable, dirTable, posTable);
    } else {
        cout << "Mode not valid...";
        exit(EXIT_FAILURE);
    }

    string locationsFileName(mode + "_locations_" + temp_comp + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
    outputLocationsFile.open(locationsFileName.c_str(), ios::out);
    string infoFileName(mode + "_info_" + temp_comp + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
    infoFile.open(infoFileName.c_str(), ios::out);

    counter = (Counters *)malloc(sizeof(Counters));
    counter->numSeeds = reads.size() * int(m / q);
    counter->numReads = reads.size();
    counter->numAcceptedSeeds = 0;

    cout << "Doing searching process..." << endl << endl;
    auto start = omp_get_wtime();

    if (searchMode.compare("all") == 0) {
        searchingReadProcess();
    } else if (searchMode.compare("exit") == 0) {
        searchingReadFoundExitProcess();
    } else {
        cout << "Invalid searching mode..." << endl;
        exit(EXIT_FAILURE);
    }

    auto end = omp_get_wtime();
    auto timeTaken = double(end - start);

    cout << "Time taken by the pigeonhole process is : " << to_string(timeTaken) << " sec" << endl << endl;
    outputLocationsFile << "Time taken by the pigeonhole process is : " << to_string(timeTaken) << " sec" << endl << endl;

    for (pair<string, vector<unsigned long long int>> readPair : forwardReadsMap) {
        vector<unsigned long long int>& temp = readPair.second;

        unordered_set<unsigned long long int> locationSet;
        for (unsigned long long int location : readPair.second)
            locationSet.insert(location);

        temp.assign(locationSet.begin(), locationSet.end());
        forwardReadsMap[readPair.first] = temp;
    }

    for (pair<string, vector<unsigned long long int>> readPair : reverseReadsMap) {
        vector<unsigned long long int>& temp = readPair.second;

        unordered_set<unsigned long long int> locationSet;
        for (unsigned long long int location : readPair.second)
            locationSet.insert(location);

        temp.assign(locationSet.begin(), locationSet.end());
        reverseReadsMap[readPair.first] = temp;
    }

//    string outFilename = "output_" + mode + "_" + temp_comp + "_" + to_string(q) + ".txt";
//    outputPossibleReadsFile.open(outFilename);
//
//    for (pair<string, vector<unsigned long long int>> readPair : readsMap) {
//        string read = readPair.first;
//
//        for (unsigned long long int location : readPair.second) {
//            if (refGenome.substr(location, m).size() == m) {
//                outputPossibleReadsFile << read << "\t" << refGenome.substr(location, m) << endl;
//            }
//        }
//    }
//
//    outputPossibleReadsFile.close();

//    counter->numAcceptedReads = forwardReadsMap.size() + reverseReadsMap.size();
//    results(forwardReadsMap, reverseReadsMap);

    cout << "Starting Bit Matrix..." << endl;
    //multiThreadedMain(forwardReadsMap, reverseReadsMap);
    multiThreadedMain();

    infoFile.close();
    outputLocationsFile.close();

    return 0;
}