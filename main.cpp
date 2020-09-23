/**
 * @file main.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the starting point of the program.
 */

#include "command.h"
#include "indexing.h"
#include "seedselector.h"
#include "bitmatrix.h"
#include "verification.h"
#include "output.h"

Genome refGenome;
vector<Read> reads;

bool isReverseAccepted;

map<unsigned long long int, vector<unsigned long long int>> minimizers;
map<long long, unsigned long long int> codeTable;
vector<unsigned long long int> dirTable;
vector<unsigned long long int> posTable;

unsigned int numSeeds;
unsigned int numReads;
unsigned int numAcceptedSeeds;
unsigned int numAcceptedReads;
unsigned int numPossibleReadLocations;
unsigned int numFilteredReadLocations;
unsigned int numVerifiedReadLocations;

string mode;
string searchMode;

string SAMFileName;

unsigned int q;
unsigned int w;
unsigned int m;
unsigned int e;
double loadFactor;

double indexRunTime;
double ssRunTime;
double bmRunTime;
double verificationRunTime;

/**
 * The starting point of the program.
 *
 * @param argc - The number of input arguments in the program.
 * @param argv - A list of input arguments in the program.
 * @return exit code
 */
int main(int argc, char *argv[]) {
    string genomeFilePath;
    string readsFilePath;
    string indexFilePath;
    string genomeFileName;

    mode = "open";
    searchMode = "exit";

    q = 20;
    e = 0;
    loadFactor = 0.9;

    isReverseAccepted = true;

    numSeeds = 0;
    numReads = 0;
    numAcceptedSeeds = 0;
    numAcceptedReads = 0;
    numPossibleReadLocations = 0;

    processingArguments(argc, argv, genomeFilePath, readsFilePath, indexFilePath, genomeFileName);

    cout << "Doing indexing process..." << endl << endl;
    auto start = omp_get_wtime();

    if (genomeFilePath.length() == 0) {
        cout << "Input reference genome was not defined..." << endl;
        exit(EXIT_FAILURE);
    }

    cout << "Reading the reference genome... " << endl << genomeFilePath << endl << endl;
    refGenome = readGenomeFile(genomeFilePath);

    if (indexFilePath.length() == 0) {
        string indexDefaultFile = mode + "_" + genomeFileName + "_" + to_string(q) + ".txt";

        ifstream indexFile(indexDefaultFile);

        if (indexFile) {
            indexFile.close();
            cout << "Index with the same file name was found by the system. ";
            readIndexFile(indexDefaultFile, minimizers, codeTable, dirTable, posTable);
        } else {
            buildIndex(genomeFileName, indexDefaultFile, minimizers, codeTable, dirTable, posTable, loadFactor);
        }
    } else {
        ifstream indexFile(indexFilePath);

        if (indexFile) {
            indexFile.close();
            readIndexFile(indexFilePath, minimizers, codeTable, dirTable, posTable);
        } else {
            cout << "Index file specified was not found by the system. ";
            buildIndex(genomeFileName, indexFilePath, minimizers, codeTable, dirTable, posTable, loadFactor);
        }
    }

    auto end = omp_get_wtime();
    auto timeTaken = double(end - start);
    cout << "Time taken by the indexing process is: " << to_string(timeTaken) << " sec" << endl << endl;
    indexRunTime = timeTaken;

    if (readsFilePath.length() == 0) {
        cout << "File for reads is not specified..." << endl;
        exit(EXIT_FAILURE);
    }

    cout << "Reading the reads... " << endl << readsFilePath << endl << endl;
    reads = readReadsFile(readsFilePath);

    w = q + q - 1;
    m = reads[0].readData.size();

    numReads = reads.size();
    numSeeds = numReads * ceil(m / (double) q);

    if (isReverseAccepted) {
        numReads *= 2;
        numSeeds *= 2;
    }

    auto filterTimeStart = omp_get_wtime();
    cout << "Doing searching process..." << endl << endl;
    start = omp_get_wtime();

    if (searchMode.compare("all") == 0) {
        seedSelectorSearchAllProcess();
    } else if (searchMode.compare("exit") == 0) {
        seedSelectorExitProcess();
    } else {
        cout << "Invalid searching mode..." << endl;
        exit(EXIT_FAILURE);
    }

    end = omp_get_wtime();
    timeTaken = double(end - start);
    ssRunTime = timeTaken;

    outputSeedSelectorResults(genomeFileName, timeTaken);
    //outputFileSeedSelectorResults(genomeFileName, timeTaken);

    /* Uncomment the line below for the pair reads output. */
    //outputPairReads(genomeFileName);

    cout << "Starting Bit Matrix..." << endl;
    bitMatrixFilterProcess();
    auto filterTimeEnd = omp_get_wtime();

    double totalFilterTime = double(filterTimeEnd - filterTimeStart);
    outputRunTimeResults(genomeFileName, indexRunTime, ssRunTime, bmRunTime, verificationRunTime, totalFilterTime);

    cout << "Starting Edlib..." << endl;
    verifyWithEdlib();
    outputPrealignmentResults();
    outputEdlibResults();

    if (SAMFileName.length() > 0) {
        outputSAMFile();
    }

    return 0;
}