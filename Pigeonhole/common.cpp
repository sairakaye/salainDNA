//
// Created by saimanalili on 28/12/2019.
//

#include "common.h"

void getDirectAddressing(string filename, vector<string> qgrams, vector<int> dirTable, vector<int> posTable) {
    ifstream directAddrFile (filename);
    string line;

    int i = 0;

    if (directAddrFile.is_open()) {
        while (getline (directAddrFile,line)) {
            if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r')) {
                line.erase(line.length()-1);
            }

            stringstream ss(line);
            string tok;

            if (i == 0) {
                while(getline(ss, tok, ' ')) {
                    if (tok.length() > 0) {
                        dirTable.push_back(stol(tok));
                    }
                }
            } else if (i == 1) {
                while(getline(ss, tok, ' ')) {
                    if (tok.length() > 0) {
                        posTable.push_back(stol(tok));
                    }
                }
            }

            i++;
        }

        directAddrFile.close();
    } else {
        cout << "File does not exist." << endl;
    }
}


void getOpenAddressing(string filename, vector<int> codeTable, vector<int> dirTable, vector<int> posTable) {
    ifstream openAddrFile (filename);
    string line;

    int i = 0;

    if (openAddrFile.is_open()) {
        while (getline (openAddrFile,line)) {
            if (!line.empty() && (line[line.length() - 1] == '\n' || line[line.length() - 1] == '\r')) {
                line.erase(line.length() - 1);
            }

            stringstream ss(line);
            string tok;

            if (i == 0) {
                while(getline(ss, tok, ' ')) {
                    if (tok.length() > 0) {
                        codeTable.push_back(stol(tok));
                    }
                }
            } else if (i == 1) {
                while(getline(ss, tok, ' ')) {
                    if (tok.length() > 0) {
                        dirTable.push_back(stol(tok));
                    }
                }
            } else if (i == 2) {
                while(getline(ss, tok, ' ')) {
                    if (tok.length() > 0) {
                        posTable.push_back(stol(tok));
                    }
                }
            }

            i++;
        }

        openAddrFile.close();
    } else {
        cout << "File does not exist." << endl;
    }
}

void generateQGrams(string prefix, vector<string>& qGrams, int k)
{
    char charactersDNA[] = {'A', 'C', 'G', 'T'};

    if (k == 0) {
        qGrams.push_back(prefix);
        return;
    }

    for (int i = 0; i < 4; i++) {
        string newPrefix;
        newPrefix = prefix + charactersDNA[i];
        generateQGrams(newPrefix, qGrams, k - 1);
    }
}

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

string readGenome(string filename) {
    string genome;

    ifstream fileGenome (filename);
    string line;

    if (fileGenome.is_open()) {
        while (getline (fileGenome,line)) {
            if (line.rfind(">", 0) == 0)
                continue;

            genome.append(line);
        }

        fileGenome.close();
    } else {
        cout << "File does not exist." << endl;
    }

    return genome;
}

vector<string> readReads(string filename) {
    vector<string> readList;

    ifstream fileRead (filename);
    string line;

    if (fileRead.is_open()) {
        while (getline (fileRead,line)) {
            if (line.rfind(">", 0) == 0)
                continue;

            if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r')) {
                line.erase(line.length()-1);
            }

            readList.push_back(line);
        }

        fileRead.close();
    } else {
        cout << "File does not exist." << endl;
        exit(EXIT_FAILURE);
    }

    return readList;
}

vector<unsigned long long> splitToInt(string str, char delimiter) {
    vector<unsigned long long> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while(getline(ss, tok, delimiter)) {
        if (tok.length() > 0) {
            internal.push_back(stol(tok));
        }
    }

    return internal;
}

map<unsigned long long, vector<unsigned long long>> getMinimizersFromFile(string filename) {
    ifstream fileMinimizer (filename);
    string line;

    map<unsigned long long, vector<unsigned long long>>  minimizers;

    if (fileMinimizer.is_open()) {
        while (getline (fileMinimizer, line)) {
            if (line.rfind(">", 0) == 0)
                continue;

            unsigned long long minimizerHashRank = stol(line.substr(0, line.find(':')));
            string locations = line.substr(line.find(':') + 1, line.length() - line.find(":"));

            minimizers[minimizerHashRank] = splitToInt(locations, ' ');
        }

        fileMinimizer.close();
    } else {
        cout << "File does not exist." << endl;
    }

    return minimizers;
}

