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

// Data of the reference genome.
Genome refGenome;

// A list of reads to be mapped by SalainDNA.
vector<Read> reads;

// It indicates if the reverse complement of each read will also be mapped or not.
bool isReverseAccepted;

// The data structure for minimizer-based indexing.
map<unsigned long long int, vector<unsigned long long int>> minimizers;
// It is where q-gram ranks are hashed into.
map<long long, unsigned long long int> codeTable;
// It contains the starting location of q-grams in the position table.
vector<unsigned long long int> dirTable;
// It is a list containing the positions of q-grams in the reference genome.
vector<unsigned long long int> posTable;

// Overall number of seeds checked by SalainDNA.
unsigned int numSeeds;
// Overall number of reads checked by SalainDNA.
unsigned int numReads;
// Overall number of accepted seeds by SalainDNA.
unsigned int numAcceptedSeeds;
// Overall number of accepted reads by SalainDNA.
unsigned int numAcceptedReads;
// Number of possible read locations found from Seed Selector.
unsigned int numPossibleReadLocations;
// Number of filtered read locations from Bit Matrix.
unsigned int numFilteredReadLocations;
// Number of verified read locations from Edlib.
unsigned int numVerifiedReadLocations;

// It indicates what mode of hash-based indexer will be used.
string mode;
// It indicates what searching mode of the seed selector will be used.
string searchMode;
// The file name of the output SAM file from SalainDNA.
string SAMFileName;

// The value of the q-gram.
unsigned int q;
// The window value to be used in minimizers Hash-based Indexer.
unsigned int w;
// The length of each read.
unsigned int m;
// The value of error threshold.
unsigned int e;
// The value of load factor used by the code table from the open addressing Hash-based Indexer.
double loadFactor;

// The runtime of the Hash-based Indexer.
double indexRunTime;
// The runtime of Seed Selector.
double ssRunTime;
// The runtime of the Bit Matrix.
double bmRunTime;
// The runtime of Edlib.
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