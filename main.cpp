#include "command.h"
#include "pigeonhole.h"
#include "bitmatrix.h"
#include "indexing.h"
#include "verification.h"

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

unsigned int numSeeds = 0;
unsigned int numReads = 0;
unsigned int numAcceptedSeeds = 0;
unsigned int numAcceptedReads = 0;
unsigned int numLocationsForward = 0;
unsigned int numLocationsReverse = 0;

ofstream infoFile;

string mode = "min";
string searchMode = "exit";

unsigned int q = 8;
unsigned int w;
unsigned int m;
unsigned int e = 0;
double loadFactor = 0.8;

int main(int argc, char *argv[]) {
    string genomeFilePath;
    string readsFilePath;
    string indexFilePath;
    string mainName;

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
            readIndexFile(indexDefaultFile, minimizers, codeTable, dirTable, posTable);
        } else {
            buildIndex(mainName, indexDefaultFile, minimizers, codeTable, dirTable, posTable, loadFactor);
        }
    } else {
        readIndexFile(indexFilePath, minimizers, codeTable, dirTable, posTable);
    }

    string infoFileName(mode + "_info_" + mainName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
    infoFile.open(infoFileName.c_str(), ios::out);

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

    cout << "Time taken by the pigeonhole process is : " << to_string(timeTaken) << " sec" << endl << endl;
    infoFile << "Time taken by the pigeonhole process is : " << to_string(timeTaken) << " sec" << endl << endl;
    removingDuplicateLocationsInEachRead();

    numAcceptedReads = forwardReadsMap.size() + reverseReadsMap.size();
    results();

    //outputPossibleReads(mainName);
    //outputPossibleLocations(mainName);

    cout << "Starting Bit Matrix..." << endl;
    multiThreadedMain();

    infoFile.close();

    return 0;
}