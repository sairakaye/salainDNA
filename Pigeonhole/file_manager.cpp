//
// Created by saimanalili on 28/12/2019.
//

#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <vector>
#include <sstream>

#include "file_manager.h"

using namespace std;

string getFileName(string filename) {
    string seperator;

    if (filename.rfind('\\') > 0) {
        seperator = "\\";
    } else if (filename.rfind('/') > 0) {
        seperator = '/';
    } else {
        return filename;
    }

    size_t sepPos = filename.rfind(seperator);

    if(sepPos != string::npos) {
        return filename.substr(sepPos + 1, filename.size() - 1);
    }

    return filename;
}

RefGenome *readGenome(string filename) {
    RefGenome *refGenome = (RefGenome *)malloc(sizeof(RefGenome));

    ifstream fileGenome (filename);
    string line;

    int genomeIndex = 0;
    int size = 0;

    refGenome->genome = (char *)malloc(1 * sizeof(char));

    if (fileGenome.is_open()) {
        while (getline (fileGenome,line)) {
            if (line.rfind(">", 0) == 0)
                continue;

            size += line.length();

            //cout << size << endl;
            //cout << line << endl;

            refGenome->genome = (char *)realloc(refGenome->genome, sizeof(char) * size);

            if (refGenome->genome == NULL) {
                cout << "Reallocation failed." << endl;
                exit(EXIT_FAILURE);
            }

            if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r')) {
                line.erase(line.length()-1);
            }

            strcpy(refGenome->genome + genomeIndex, &line[0]);
            genomeIndex = size;
        }

        refGenome->genome = (char *)realloc(refGenome->genome, sizeof(char) * size + 1);

        if (refGenome->genome == NULL) {
            cout << "Reallocation failed." << endl;
            exit(EXIT_FAILURE);
        } else {
            *(refGenome->genome + size) = '\0';
        }

        refGenome->length = size;

        fileGenome.close();
    } else {
        cout << "File does not exist." << endl;
    }

    return refGenome;
}

ReadList *readReads(string filename) {
    ReadList *readList = (ReadList *)malloc(sizeof(ReadList));

    char **reads = (char **)malloc(sizeof(char *));

    ifstream fileRead (filename);
    string line;

    int num = 1;

    if (fileRead.is_open()) {
        int size = 0;

        while (getline (fileRead,line)) {
            if (line.rfind(">", 0) == 0)
                continue;

            reads = (char **)realloc(reads, num * sizeof(char *));

            if (reads == NULL) {
                cout << "Reallocation failed." << endl;
                exit(EXIT_FAILURE);
            }

            size = line.length();

            *(reads + (num - 1)) = (char *)malloc((size + 1) * sizeof(char));

            if (*(reads + (num - 1)) == NULL) {
                cout << "Memory allocation failed." << endl;
                exit(EXIT_FAILURE);
            }

            if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r')) {
                line.erase(line.length()-1);
            }

            strcpy(*(reads + (num - 1)), &line[0]);

            num++;
        }

        fileRead.close();
    } else {
        cout << "File does not exist." << endl;
        return NULL;
    }

    readList->reads = reads;
    readList->size = num - 1;

    return readList;
}

vector<int> splitToInt(string str, char delimiter) {
    vector<int> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while(getline(ss, tok, delimiter)) {
        if (tok.length() > 0) {
            internal.push_back(stoi(tok));
        }
    }

    return internal;
}

map<int, vector<int>> getMinimizersFromFile(string filename) {
    ifstream fileMinimizer (filename);
    string line;

    map<int, vector<int>>  minimizers;

    if (fileMinimizer.is_open()) {
        while (getline (fileMinimizer, line)) {
            if (line.rfind(">", 0) == 0)
                continue;

            int minimizerHashRank = stoi(line.substr(0, line.find(':')));
            string locations = line.substr(line.find(':') + 1, line.length() - line.find(":"));

            minimizers[minimizerHashRank] = splitToInt(locations, ' ');
        }

        fileMinimizer.close();
    } else {
        cout << "File does not exist." << endl;
    }

    return minimizers;
}

