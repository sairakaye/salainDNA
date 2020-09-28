/**
 * @file command.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the implementation of the processing of the arugments.
 */

#include "command.h"

/**
 * It extracts the file name from the given file path.
 *
 * @param filePath - The file path where the file name will be extracted.
 * @return the file name
 */
string getFileName(string filePath) {
    char separator = '/';

    #ifdef _WIN32
        separator = '\\';
    #endif

    unsigned int i = filePath.rfind(separator, filePath.length());

    if (i != string::npos) {
        return (filePath.substr(i+1, filePath.length() - i));
    }

    return ("");
}

/**
 * It prints the general details of SalainDNA. The details contain the name of the program and its features.
 *
 */
void generalDetails() {
    cout << "SalainDNA: A Pre-alignment Filter for DNA Read Mapping on a Multi-core Environment" << endl << endl;

    cout << "The program contains:" << endl;
    cout << "\t - Indexing (Direct Addressing, Open Addressing, Minimizers)" << endl;
    cout << "\t - Seed Selector using Generalized Pigeonhole Principle" << endl;
    cout << "\t - Bit Matrix for further filteration" << endl;
    cout << "\t - Edlib for Verification" << endl;
}

/**
 * It prints the commands or arguments that the user can use to run the program.
 *
 */
void helpCommands() {
    cout << "Arguments:" << endl;
    cout << "\t -h - view the commands" << endl;
    cout << "\t -v - view the version" << endl;
    cout << "\t -q - set the q-gram for the seed and the indexing (default is 20)" << endl;
    cout << "\t -g - file path of the genome file for generating index" << endl;
    cout << "\t -ir - file path of the input reads for mapping" << endl;
    cout << "\t -i - file path for the index file" << endl;
    cout << "\t -m - mode of indexing (dir - direct addressing, open - open addressing, min - minimizers) (default is open)" << endl;
    cout << "\t -l - load factor for open addressing indexing mode (default is 0.9)" << endl;
    cout << "\t -e - error threshold (default is 0)" << endl;
    cout << "\t -s - search mode for seed selector (all - find locations in each seed, exit - when a location is found in the seed, immediately exit)" << endl;
    cout << "\t -rev - find using reverse complement (0 - false, 1 - true) (default is 1)" << endl;
    cout << "\t -out - output file for the SAM file" << endl;
}

/**
 * It processes the command or arguments input in the program by the user.
 *
 * @param argc - The number of input arguments in the program.
 * @param argv - A list of input arguments in the program.
 * @param genomeFilePath - The file path of the reference genome.
 * @param readsFilePath - The file path of the reads to be mapped by SalainDNA.
 * @param indexFilePath - The file path of the index.
 * @param genomeFileName - The filename of the reference genome.
 */
void processingArguments(int argc, char *argv[], string &genomeFilePath, string &readsFilePath, string &indexFilePath, string &genomeFileName) {
    if (argc == 1) {
        generalDetails();
        cout << endl;
        helpCommands();
        exit(EXIT_SUCCESS);
    } else {
        for (int i = 1; i < argc; ++i) {
            if (string(argv[i]) == "-h") {
                helpCommands();
                exit(EXIT_SUCCESS);
            } else if (string(argv[i]) == "-v") {
                cout << "Version 1.0.0" << endl;
                exit(EXIT_SUCCESS);
            } else if (string(argv[i]) == "-q") {
                try {
                    q = stoi(argv[i + 1]);
                } catch (exception &err) {
                    cout << err.what() << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (string(argv[i]) == "-g") {
                genomeFilePath = string(argv[i + 1]);
                genomeFileName = getFileName(genomeFilePath);
                genomeFileName = genomeFileName.substr(0, genomeFileName.find_last_of("."));
            } else if (string(argv[i]) == "-ir") {
                readsFilePath = string(argv[i + 1]);
            } else if (string(argv[i]) == "-i") {
                indexFilePath = string(argv[i + 1]);
            } else if (string(argv[i]) == "-m") {
                mode = string(argv[i + 1]);
            } else if (string(argv[i]) == "-e") {
                try {
                    e = stoi(argv[i + 1]);
                } catch (exception &err) {
                    cout << err.what() << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (string(argv[i]) == "-s") {
                searchMode = string(argv[i + 1]);
            } else if (string(argv[i]) == "-l") {
                try {
                    loadFactor = stod(argv[i + 1]);
                } catch (exception &err) {
                    cout << err.what() << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (string(argv[i]) == "-rev") {
                int value;

                try {
                    value = stoi(argv[i+1]);

                    if (value == 0) {
                        isReverseAccepted = false;
                    } else if (value == 1) {
                        isReverseAccepted = true;
                    } else {
                        cout << "Input either 1 for true or 0 for false." << endl;
                        exit(EXIT_FAILURE);
                    }
                } catch (exception &err) {
                    cout << "Invalid input. Input either 1 for true or 0 for false." << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (string(argv[i]) == "-out") {
                SAMFileName = string(argv[i + 1]);
            }
        }
    }
}