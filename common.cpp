//
// Created by saimanalili on 25/02/2020.
//

#include "common.h"

vector<pair<string, int>> alphabetRef = { {"A", 0}, {"C",1}, {"G",2}, {"T", 3} };

unsigned long long int extractRanking(string kMer) {
    string binary;
    int rankValue;

    for (int i = 0; i<kMer.length(); i++){
        for (int j = 0; j< alphabetRef.size(); j++){
            if (kMer.at(i) + string() == alphabetRef.at(j).first){
                rankValue = alphabetRef.at(j).second;
                binary.append(bitset<2>(rankValue).to_string());
            }
        }
    }

    unsigned long long int decimal = strtoull(binary.c_str(), nullptr, 2);
    return decimal;
}

uint64_t inthash_64(uint64_t key, uint64_t mask) {
    key = (~key + (key << 21)) & mask;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

Genome readGenomeFile(string filename) {
    Genome genome;
    string genomeData;
    string genomeName;
    ifstream fileGenome (filename);
    string line;

    if (fileGenome.is_open()) {
        while (getline (fileGenome,line)) {
            if (line.rfind(">", 0) == 0) {
                genomeName = line.substr(1, line.find(' ') - 1);
                continue;
            }

            if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r')) {
                line.erase(line.length()-1);
            }

            genomeData.append(line);
        }

        for (int j = 0; j < genomeData.length(); j++) {
            if (genomeData.at(j) != 'A' && genomeData.at(j) != 'C' && genomeData.at(j) != 'G' && genomeData.at(j) != 'T')
                replace(genomeData.begin(), genomeData.end(), genomeData.at(j), '\0');
        }

        genome.genomeName = genomeName;
        genome.genomeData = genomeData;

        fileGenome.close();
    } else {
        cout << "File does not exist." << endl;
        cout << "File name: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    return genome;
}

vector<Read> readReadsFile(string filename) {
    vector<Read> readList;
    ifstream fileRead(filename);
    string line;
    string read;
    string readName;

    if (fileRead.is_open()) {
        while (getline (fileRead,line)) {
            if (line.rfind(">", 0) == 0) {
                if (read.size() > 0 && readName.size() > 0) {
                    Read readStruct;
                    readStruct.readName = readName;
                    readStruct.readData = read;
                    //readMap[read] = readStruct;
                    readList.push_back(readStruct);
                    read = "";
                    //numReads++;
                }

                readName = line.substr(1, line.find(' ') - 1);
                continue;
            }

            if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r')) {
                line.erase(line.length()-1);
            }

            read.append(line);

            /*
            readList.push_back(line);
            readsLabelMap[line] = readName;
            */
        }

        if (read.size() > 0 && readName.size() > 0) {
            Read readStruct;
            readStruct.readName = readName;
            readStruct.readData = read;
            //readMap[read] = readStruct;
            readList.push_back(readStruct);
        }

        fileRead.close();
    } else {
        cout << "File does not exist." << endl;
        cout << "File name: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    return readList;
}

string reverseComplement(string read) {
    string reverseRead = read;
    reverse(reverseRead.begin(), reverseRead.end());

    for (int i = 0; i < reverseRead.length(); i++) {
        switch (reverseRead[i]) {
            case 'A':
                reverseRead[i] = 'T';
                break;
            case 'C':
                reverseRead[i] = 'G';
                break;
            case 'G':
                reverseRead[i] = 'C';
                break;
            case 'T':
                reverseRead[i] = 'A';
                break;
        }
    }

    return reverseRead;
}

void getDirectAddressing(string filename, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable) {
    ifstream directAddrFile(filename);
    string line;

    bool isGettingDir = false;
    bool isGettingPos = false;

    if (directAddrFile.is_open()) {
        while (getline (directAddrFile,line)) {
            if (!line.empty() && (line[line.length()-1] == '\n' || line[line.length()-1] == '\r')) {
                line.erase(line.length()-1);
            }

            stringstream ss(line);
            string tok;

            if (line.compare("dir") == 0) {
                isGettingDir = true;
                isGettingPos = false;
                continue;
            } else if (line.compare("pos") == 0) {
                isGettingDir = false;
                isGettingPos = true;
                continue;
            } else {
                if (isGettingDir) {
                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0) {
                            dirTable.push_back(strtoull(tok.c_str(), nullptr, 10));
                        }
                    }
                } else if (isGettingPos) {
                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0) {
                            posTable.push_back(strtoull(tok.c_str(), nullptr, 10));
                        }
                    }
                }
            }
        }

        directAddrFile.close();
    } else {
        cout << "File does not exist." << endl;
        cout << "File name: " << filename << endl;
        exit(EXIT_FAILURE);
    }
}

void getOpenAddressing(string filename, map<long long, unsigned long long int>& codeTable, vector<unsigned long long int>& dirTable, vector<unsigned long long int>& posTable) {
    ifstream openAddrFile(filename);
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
                    long long codeRank;
                    unsigned long long int codeIndex;
                    int i = 0;

                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0 && tok.compare("") != 0) {
                            if (i == 0) {
                                codeRank = strtoll(tok.c_str(), nullptr, 10);

                                if (codeRank == -1) {
                                    break;
                                }
                            } else if (i == 1) {
                                codeIndex = strtoull(tok.c_str(), nullptr, 10);
                            }
                        }

                        i++;
                    }

                    if (codeRank != -1) {
                        codeTable.insert(pair<long long, unsigned long long int>(codeRank, codeIndex));
                    }
                } else if (isGettingDir) {
                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0) {
                            dirTable.push_back(strtoull(tok.c_str(), nullptr, 10));
                        }
                    }
                } else if (isGettingPos) {
                    while(getline(ss, tok, ' ')) {
                        if (tok.length() > 0) {
                            posTable.push_back(strtoull(tok.c_str(), nullptr, 10));
                        }
                    }
                }
            }
        }

        openAddrFile.close();
    } else {
        cout << "File does not exist." << endl;
        cout << "File name: " << filename << endl;
        exit(EXIT_FAILURE);
    }
}

vector<unsigned long long int> splitToInt(string str, char delimiter) {
    vector<unsigned long long int> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while(getline(ss, tok, delimiter)) {
        if (tok.length() > 0) {
            internal.push_back(strtoull(tok.c_str(), nullptr, 10));
        }
    }

    return internal;
}

map<unsigned long long int, vector<unsigned long long int>> getMinimizers(string filename) {
    ifstream fileMinimizer(filename);
    string line;

    map<unsigned long long int, vector<unsigned long long int>> minimizers;

    if (fileMinimizer.is_open()) {
        while (getline (fileMinimizer, line)) {
            if (line.rfind(">", 0) == 0)
                continue;

            unsigned long long int minimizerHashRank = strtoull(line.substr(0, line.find(':')).c_str(), nullptr, 10);
            string locations = line.substr(line.find(':') + 1, line.length() - line.find(":"));

            minimizers[minimizerHashRank] = splitToInt(locations, ' ');
        }

        fileMinimizer.close();
    } else {
        cout << "File does not exist." << endl;
        cout << "File name: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    return minimizers;
}

/*
void processingPossibleReadsForBitmatrix() {
    for (pair<string, vector<unsigned long long int>> readPair : possibleReadsMap) {
        vector<unsigned long long int> &temp = readPair.second;

        unordered_set<unsigned long long int> locationSet;
        for (unsigned long long int location : readPair.second) {
            if (refGenome.substr(location, m).size() == m) {
                locationSet.insert(location);
            }
        }

        temp.assign(locationSet.begin(), locationSet.end());
        possibleReadsMap[readPair.first] = temp;

        numLocations += temp.size();
    }
}
*/