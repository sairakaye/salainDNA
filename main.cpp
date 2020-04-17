#include "command.h"
#include "pigeonhole.h"
#include "bitmatrix.h"
#include "indexing.h"
#include "output.h"
//#include "verification.h"

string refGenome;
map<string, string> reads;

map<string, vector<unsigned long long int>> forwardReadsMap;
map<string, vector<unsigned long long int>> reverseReadsMap;
map<string, vector<unsigned long long int>> filteredReadsMap;

map<unsigned long long int, vector<unsigned long long int>> minimizers;
map<long long, unsigned long long int> codeTable;
vector<unsigned long long int> dirTable;
vector<unsigned long long int> posTable;

unsigned int numSeeds;
unsigned int numReads;
unsigned int numAcceptedSeeds;
unsigned int numAcceptedReads;
unsigned int numLocationsForward;
unsigned int numLocationsReverse;

string mode;
string searchMode;

unsigned int q;
unsigned int w;
unsigned int m;
unsigned int e;
double loadFactor = 0.8;

int main(int argc, char *argv[]) {
    string genomeFilePath;
    string readsFilePath;
    string indexFilePath;
    string mainName;

    mode = "min";
    searchMode = "exit";

    q = 8;
    e = 0;

    numSeeds = 0;
    numReads = 0;
    numAcceptedSeeds = 0;
    numAcceptedReads = 0;
    numLocationsForward = 0;
    numLocationsReverse = 0;

    processingArguments(argc, argv, genomeFilePath, readsFilePath, indexFilePath, mainName);

    cout << "Reading the reference genome... " << endl << genomeFilePath << endl << endl;
    refGenome = readGenomeFile(genomeFilePath);
    cout << "Reading the reads... " << endl << readsFilePath << endl << endl;
    reads = readReadsFile(readsFilePath);

    w = q + q - 1;
    m = reads.begin()->first.size();

    if (indexFilePath.length() == 0) {
        string indexDefaultFile = mode + "_" + mainName + "_" + to_string(q) + ".txt";

        ifstream indexFile(indexDefaultFile);

        if (indexFile) {
            indexFile.close();
            readIndexFile(indexDefaultFile, minimizers, codeTable, dirTable, posTable);
        } else {
            buildIndex(mainName, indexDefaultFile, minimizers, codeTable, dirTable, posTable, loadFactor);
        }
    } else {
        readIndexFile(indexFilePath, minimizers, codeTable, dirTable, posTable);
    }

    numSeeds = reads.size() * int(m / q);
    numReads = reads.size();

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

    removingDuplicateLocationsInEachRead();
    numAcceptedReads = forwardReadsMap.size() + reverseReadsMap.size();

    outputSeedSelectorResults(mainName, timeTaken);

    //outputPossibleReads(mainName);
    //outputPossibleLocations(mainName);

    cout << "Starting Bit Matrix..." << endl;
    multiThreadedMain();

    return 0;
}