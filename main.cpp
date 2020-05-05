#include "command.h"
#include "pigeonhole.h"
#include "bitmatrix.h"
#include "indexing.h"
#include "output.h"
//#include "verification.h"

string refGenome;
vector<string> reads;
map<string, string> readsLabelMap;

//map<string, vector<unsigned long long int>> forwardReadsMap;
//map<string, vector<unsigned long long int>> reverseReadsMap;
map<string, vector<unsigned long long int>> possibleReadsMap;
map<string, vector<unsigned long long int>> filteredReadsMap;

//vector<PossibleRead> possibleReadsVector;
//vector<PossibleRead> filteredReadsVector;

bool isReverseAccepted;

map<unsigned long long int, vector<unsigned long long int>> minimizers;
map<long long, unsigned long long int> codeTable;
vector<unsigned long long int> dirTable;
vector<unsigned long long int> posTable;

unsigned int numSeeds;
unsigned int numReads;
unsigned int numAcceptedSeeds;
unsigned int numAcceptedReads;
unsigned int numLocations;

string mode;
string searchMode;

unsigned int q;
unsigned int w;
unsigned int m;
unsigned int e;
double loadFactor;

int main(int argc, char *argv[]) {
    string genomeFilePath;
    string readsFilePath;
    string indexFilePath;
    string mainName;

    mode = "min";
    searchMode = "exit";

    q = 12;
    e = 0;
    loadFactor = 0.9;

    isReverseAccepted = false;

    numSeeds = 0;
    numReads = 0;
    numAcceptedSeeds = 0;
    numAcceptedReads = 0;
    numLocations = 0;

    processingArguments(argc, argv, genomeFilePath, readsFilePath, indexFilePath, mainName);

    if (genomeFilePath.length() == 0) {
        cout << "Input reference genome was not defined..." << endl;
        exit(EXIT_FAILURE);
    }

    cout << "Reading the reference genome... " << endl << genomeFilePath << endl << endl;
    refGenome = readGenomeFile(genomeFilePath);

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
        ifstream indexFile(indexFilePath);

        if (indexFile) {
            indexFile.close();
            smatch m;
            regex indexRegex("([min|dir|open]+)_([a-zA-z0-9_]*)_([0-9]*).txt");

            if (regex_search(indexFilePath, m, indexRegex)) {
                if (mode == m.str(1) && to_string(q) == m.str(3)) {
                    mainName = m.str(2);
                    readIndexFile(indexFilePath, minimizers, codeTable, dirTable, posTable);
                } else {
                    cout << "File does not match mode and q set..." << endl;
                    exit(EXIT_FAILURE);
                }
            } else {
                cout << "Index file name did not match system defined file name..." << endl;
                exit(EXIT_FAILURE);
            }
        } else {
            cout << "File does not exist." << endl;
            cout << "File name: " << indexFilePath << endl << endl;
            buildIndex(mainName, indexFilePath, minimizers, codeTable, dirTable, posTable, loadFactor);
        }
    }

    if (readsFilePath.length() == 0) {
        cout << "Indexing done..." << endl;
        exit(EXIT_FAILURE);
    }

    cout << "Reading the reads... " << endl << readsFilePath << endl << endl;
    reads = readReadsFile(readsFilePath);

    w = q + q - 1;
    m = reads[0].size();

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

    numAcceptedReads = possibleReadsMap.size();

    outputSeedSelectorResults(mainName, timeTaken);

    //outputPossibleReads(mainName);
    //outputPossibleLocations(mainName);

    cout << "Starting Bit Matrix..." << endl;
    multiThreadedMain();

    return 0;
}