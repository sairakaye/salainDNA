//
// Created by saimanalili on 17/04/2020.
//

#include "output.h"

void outputPossibleReads(string& mainName) {
    ofstream outputPossibleReadsFile;

    string outFilename = "output_possible_reads_" + mode + "_" + mainName + "_" + to_string(q) + ".txt";
    outputPossibleReadsFile.open(outFilename);

    for (pair<string, vector<unsigned long long int>> readPair : possibleReadsMap) {
        string read = readPair.first;

        for (unsigned long long int location : readPair.second) {
            if (refGenome.substr(location, m).size() == m) {
                outputPossibleReadsFile << read << "\t" << refGenome.substr(location, m) << endl;
            }
        }
    }

    /*

    for (pair<string, vector<unsigned long long int>> readPair : forwardReadsMap) {
        string read = readPair.first;

        for (unsigned long long int location : readPair.second) {
            if (refGenome.substr(location, m).size() == m) {
                outputPossibleReadsFile << read << "\t" << refGenome.substr(location, m) << endl;
            }
        }
    }

    for (pair<string, vector<unsigned long long int>> readPair : reverseReadsMap) {
        string read = readPair.first;

        for (unsigned long long int location : readPair.second) {
            if (refGenome.substr(location, m).size() == m) {
                outputPossibleReadsFile << read << "\t" << refGenome.substr(location, m) << endl;
            }
        }
    }
    */

    outputPossibleReadsFile.close();
}

void outputPossibleLocations(string& mainName) {
    ofstream outputPossibleLocations;

    string outFilename = "output_possible_locations_" + mode + "_" + mainName + "_" + to_string(q) + ".txt";
    outputPossibleLocations.open(outFilename);

    for (pair<string, vector<unsigned long long int>> readPair : possibleReadsMap) {
        string read = readPair.first;

        outputPossibleLocations << read << endl;

        for (unsigned long long int location : readPair.second) {
            if (refGenome.substr(location, m).size() == m) {
                outputPossibleLocations << location << endl;
            }
        }

        outputPossibleLocations << endl;
    }

    /*
    for (pair<string, vector<unsigned long long int>> readPair : forwardReadsMap) {
        string read = readPair.first;

        outputPossibleLocations << read << endl;

        for (unsigned long long int location : readPair.second) {
            if (refGenome.substr(location, m).size() == m) {
                outputPossibleLocations << location << endl;
            }
        }

        outputPossibleLocations << endl;
    }

    for (pair<string, vector<unsigned long long int>> readPair : reverseReadsMap) {
        string read = readPair.first;

        outputPossibleLocations << read << endl;

        for (unsigned long long int location : readPair.second) {
            if (refGenome.substr(location, m).size() == m) {
                outputPossibleLocations << location << endl;
            }
        }

        outputPossibleLocations << endl;
    }
    */

    outputPossibleLocations.close();
}

void outputSeedSelectorResults(string& mainName, double timeTaken) {
    //ofstream infoFile;
    //string infoFileName(mode + "_info_" + mainName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
    //infoFile.open(infoFileName.c_str(), ios::out);

    cout << "Time taken by the pigeonhole process is : " << to_string(timeTaken) << " sec" << endl << endl;
    //infoFile << "Time taken by the pigeonhole process is : " << to_string(timeTaken) << " sec" << endl << endl;

    cout << "Number of seeds checked: " << to_string(numSeeds) << endl;
    cout << "Number of reads checked: " << to_string(numReads) << endl << endl;
    //infoFile << "Number of seeds checked: " << to_string(numSeeds) << endl;
    //infoFile << "Number of reads checked: " << to_string(numReads) << endl << endl;

    cout << "Number of accepted seeds: " << to_string(numAcceptedSeeds) << endl;
    cout << "Number of accepted reads: " << to_string(numAcceptedReads) << endl << endl;
    //infoFile << "Number of accepted seeds: " << to_string(numAcceptedSeeds) << endl;
    //infoFile << "Number of accepted reads: " << to_string(numAcceptedReads) << endl << endl;


    cout << "Number of possible read locations found: " + to_string(numLocations) << endl << endl;
    //infoFile << "Number of possible read locations found from forward: " + to_string(numLocations) << endl;

    //infoFile.close();
}

void outputPrealignmentResults() {

}