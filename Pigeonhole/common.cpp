//
// Created by saimanalili on 28/12/2019.
//

#include "common.h"

void getDirectAddressing(string filename, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable) {
    ifstream directAddrFile (filename);
    string line;

    int isGettingDir = 0;
    int isGettingPos = 0;

    if (directAddrFile.is_open()) {
        while (getline (directAddrFile,line)) {
            if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r')) {
                line.erase(line.length()-1);
            }

            stringstream ss(line);
            string tok;

            if (line.compare("dir") == 0) {
                isGettingDir = 1;
                isGettingPos = 0;
                continue;
            } else if (line.compare("pos") == 0) {
                isGettingDir = 0;
                isGettingPos = 1;
                continue;
            } else {
                if (isGettingDir == 1) {
                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0) {
                            dirTable.push_back(stol(tok));
                        }
                    }
                } else if (isGettingPos == 1) {
                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0) {
                            posTable.push_back(stol(tok));
                        }
                    }
                }
            }
        }

        directAddrFile.close();
    } else {
        cout << "File does not exist." << endl;
    }
}


void getOpenAddressing(string filename, map<long long, unsigned long long int>& codeTable, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable) {
    ifstream openAddrFile (filename);
    string line;

    bool isGettingCode = false;
    bool isGettingDir = false;
    bool isGettingPos = false;

    if (openAddrFile.is_open()) {
        while (getline (openAddrFile,line)) {
            if (!line.empty() && (line[line.length() - 1] == '\n' || line[line.length() - 1] == '\r')) {
                line.erase(line.length() - 1);
            }

            stringstream ss(line);
            string tok;
            if (line.compare("code") == 0) {
                isGettingCode = true;
                isGettingDir = false;
                isGettingPos = false;
                continue;
            } else if (line.compare("dir") == 0) {
                isGettingCode = false;
                isGettingDir = true;
                isGettingPos = false;
                continue;
            } else if (line.compare("pos") == 0) {
                isGettingCode = false;
                isGettingDir = false;
                isGettingPos = true;
                continue;
            } else {
                if (isGettingCode) {
                    long long getThis;
                    unsigned long long int num;
                    int i = 0;

                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0 && tok.compare("") != 0) {
                            if (i == 0) {
                                getThis = stoll(tok);

                                if (getThis == -1) {
                                    break;
                                }
                            } else if (i == 1) {
                                num = stoll(tok);
                            }
                        }

                        i++;
                    }

                    if (getThis != -1) {
                        codeTable.insert(pair<long long, unsigned long long int>(getThis, num));
                    }
                } else if (isGettingDir) {
                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0) {
                            dirTable.push_back(stoll(tok));
                        }
                    }
                } else if (isGettingPos) {
                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0) {
                            posTable.push_back(stoll(tok));
                        }
                    }
                }
            }
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
        exit(EXIT_FAILURE);
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

vector<unsigned long long int> splitToInt(string str, char delimiter) {
    vector<unsigned long long int> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while(getline(ss, tok, delimiter)) {
        if (tok.length() > 0) {
            internal.push_back(stol(tok));
        }
    }

    return internal;
}

map<unsigned long long int, vector<unsigned long long int>> getMinimizersFromFile(string filename) {
    ifstream fileMinimizer (filename);
    string line;

    map<unsigned long long int, vector<unsigned long long int>> minimizers;

    if (fileMinimizer.is_open()) {
        while (getline (fileMinimizer, line)) {
            if (line.rfind(">", 0) == 0)
                continue;

            unsigned long long int minimizerHashRank = stol(line.substr(0, line.find(':')));
            string locations = line.substr(line.find(':') + 1, line.length() - line.find(":"));

            minimizers[minimizerHashRank] = splitToInt(locations, ' ');
        }

        fileMinimizer.close();
    } else {
        cout << "File does not exist." << endl;
    }

    return minimizers;
}
