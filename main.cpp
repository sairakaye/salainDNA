#include "command.h"
#include "pigeonhole.h"
#include "bitmatrix.h"
#include "indexing.h"
#include "verification.h"

using namespace std::chrono;

string refGenome;
vector<string> reads;
map<string, vector<string>> readsLabelMap;

map<string, vector<unsigned long long int>> forwardReadsMap;
map<string, vector<unsigned long long int>> reverseReadsMap;
map<string, vector<unsigned long long int>> filteredReadsMap;

map<unsigned long long int, vector<unsigned long long int>> minimizers;
map<long long, unsigned long long int> codeTable;
vector<unsigned long long int> dirTable;
vector<unsigned long long int> posTable;

ofstream outputPossibleReadsFile;
ofstream outputLocationsFile;
ofstream infoFile;

// Change. Kawawa memory haha.
Counters *counter;

string mode;
string searchMode;
unsigned int q;
unsigned int w;
unsigned int m;
unsigned int e;

int main(int argc, char *argv[]) {
    string genomeFilePath;
    string readsFilePath;
    string indexFilePath;
    string mainName;

    q = 8;
    e = 0;
    mode = "min";
    searchMode = "exit";

    processingArguments(argc, argv, genomeFilePath, readsFilePath, indexFilePath, mainName);

    cout << "Reading the reference genome... " << endl << genomeFilePath << endl << endl;
    refGenome = readGenomeFile(genomeFilePath);
    cout << "Reading the reads... " << endl << readsFilePath << endl << endl;
    reads = readReadsFile(readsFilePath);

    w = q + q - 1;
    m = reads[0].size();

    if (indexFilePath.length() == 0) {
        string indexDefaultFile = mode + "_" + mainName + "_" + to_string(q) + ".txt";

        ifstream indexFile(indexDefaultFile);

        if (indexFile) {
            readIndexFile(mode, indexDefaultFile, minimizers, codeTable, dirTable, posTable);
        } else {
            buildIndex(mode, mainName, indexDefaultFile, minimizers, codeTable, dirTable, posTable);
        }
    } else {
        readIndexFile(mode, indexFilePath, minimizers, codeTable, dirTable, posTable);
    }

    string locationsFileName(mode + "_locations_" + mainName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
    outputLocationsFile.open(locationsFileName.c_str(), ios::out);
    string infoFileName(mode + "_info_" + mainName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
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

    results(forwardReadsMap, reverseReadsMap);

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