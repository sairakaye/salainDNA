//
// Created by saimanalili on 28/12/2019.
//

#include "common.h"

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

