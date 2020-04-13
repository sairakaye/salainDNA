//
// Created by Sai Manalili on 4/8/2020.
//

#include "command.h"

string getFileName(string filePath) {

    char sep = '/';

#ifdef _WIN32
    sep = '\\';
#endif

    unsigned int i = filePath.rfind(sep, filePath.length());

    if (i != string::npos) {
        return (filePath.substr(i+1, filePath.length() - i));
    }

    return ("");
}

void processingArguments(int argc, char *argv[], string &genomeFilePath, string &readsFilePath, string &indexFilePath, string &mainName) {
    if (argc == 1) {
        cout << "Print the details of the program." << endl;
        exit(EXIT_SUCCESS);
    } else {
        for (int i = 1; i < argc; ++i) {
            if (string(argv[i]) == "-h") {
                cout << "Help commands" << endl;
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
                mainName = getFileName(genomeFilePath);
                mainName = mainName.substr(0, mainName.find_last_of("."));
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
            }
        }
    }
}
