//
// Created by Sai Manalili on 4/8/2020.
//

#include "command.h"

// Tamang try catch lang hehe mamaya.
void processingArguments(int argc, char *argv[], string &genomeFileName, string &readsFilename, string &indexFilename, string &mainName) {
    if (argc == 1) {
        cout << "Print the deets." << endl;
        exit(EXIT_SUCCESS);
    } else {
        for (int i = 1; i < argc; ++i) {
            if (string(argv[i]) == "-h") {
                cout << "Help deets" << endl;
                exit(EXIT_SUCCESS);
            } else if (string(argv[i]) == "-v") {
                cout << "Version deets" << endl;
                exit(EXIT_SUCCESS);
            } else if (string(argv[i]) == "-q") {
                try {
                    q = stoi(argv[i + 1]);
                } catch (exception &err) {
                    cout << err.what() << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (string(argv[i]) == "-g") {
                genomeFileName = string(argv[i + 1]);
            } else if (string(argv[i]) == "-ir") {
                readsFilename = string(argv[i + 1]);
            } else if (string(argv[i]) == "-i") {
                indexFilename = string(argv[i + 1]);
            } else if (string(argv[i]) == "-m") {
                mode = string(argv[i + 1]);
            } else if (string(argv[i]) == "-temp") {
                temp_comp = string(argv[i + 1]);
            } else if (string(argv[i]) == "-e") {
                try {
                    e = stoi(argv[i + 1]);
                } catch (exception &err) {
                    cout << err.what() << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (string(argv[i]) == "-s") {
                searchMode = string(argv[i + 1]);
            }
        }
    }
}
