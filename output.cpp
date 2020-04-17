//
// Created by saimanalili on 17/04/2020.
//

#include "output.h"

void outputPossibleReads(string& mainName) {
    ofstream outputPossibleReadsFile;

    string outFilename = "output_possible_reads_" + mode + "_" + mainName + "_" + to_string(q) + ".txt";
    outputPossibleReadsFile.open(outFilename);

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

    outputPossibleReadsFile.close();
}

void outputPossibleLocations(string& mainName) {
    ofstream outputPossibleLocations;

    string outFilename = "output_possible_locations_" + mode + "_" + mainName + "_" + to_string(q) + ".txt";
    outputPossibleLocations.open(outFilename);

    for (pair<string, vector<unsigned long long int>> readPair : forwardReadsMap) {
        string read = readPair.first;

        outputPossibleLocations << read + " - ";

        for (unsigned long long int location : readPair.second) {
            if (refGenome.substr(location, m).size() == m) {
                outputPossibleLocations << location + " ";
            }
        }

        outputPossibleLocations << endl;
    }

    for (pair<string, vector<unsigned long long int>> readPair : reverseReadsMap) {
        string read = readPair.first;

        outputPossibleLocations << read + " - ";

        for (unsigned long long int location : readPair.second) {
            if (refGenome.substr(location, m).size() == m) {
                outputPossibleLocations << read << "\t" << refGenome.substr(location, m) << endl;
            }
        }

        outputPossibleLocations << endl;
    }

    outputPossibleLocations.close();
}

void outputSeedSelectorResults(string& mainName, double timeTaken) {
    ofstream infoFile;
    string infoFileName(mode + "_info_" + mainName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
    infoFile.open(infoFileName.c_str(), ios::out);

    cout << "Time taken by the pigeonhole process is : " << to_string(timeTaken) << " sec" << endl << endl;
    infoFile << "Time taken by the pigeonhole process is : " << to_string(timeTaken) << " sec" << endl << endl;

    cout << "Number of seeds checked: " << to_string(numSeeds) << endl;
    cout << "Number of reads checked: " << to_string(numReads) << endl << endl;
    infoFile << "Number of seeds checked: " << to_string(numSeeds) << endl;
    infoFile << "Number of reads checked: " << to_string(numReads) << endl << endl;

    cout << "Number of accepted seeds: " << to_string(numAcceptedSeeds) << endl;
    cout << "Number of accepted reads: " << to_string(numAcceptedReads) << endl << endl;
    infoFile << "Number of accepted seeds: " << to_string(numAcceptedSeeds) << endl;
    infoFile << "Number of accepted reads: " << to_string(numAcceptedReads) << endl << endl;

//    vector<unsigned long long int> forwardFound;
//    vector<unsigned long long int> reverseFound;
//
//    cout << "Putting forward and reverse locations for counting..." << endl;
//    for (pair<string, vector<unsigned long long int>> readPair : forwardReadsMap) {
//        for (unsigned long long int location : readPair.second) {
//            forwardFound.push_back(location);
//        }
//    }
//
//    for (pair<string, vector<unsigned long long int>> readPair : reverseReadsMap) {
//        for (unsigned long long int location : readPair.second) {
//            reverseFound.push_back(location);
//        }
//    }
//
//    cout << "Removing duplicate locations mapped..." << endl;
//    unordered_set<unsigned long long int> duplicateRemoverSet;
//
//    for (unsigned long long int location : forwardFound) {
//        duplicateRemoverSet.insert(location);
//    }
//
//    forwardFound.assign(duplicateRemoverSet.begin(), duplicateRemoverSet.end());
//    duplicateRemoverSet.clear();
//
//    for (unsigned long long int location : reverseFound) {
//        duplicateRemoverSet.insert(location);
//    }
//
//    reverseFound.assign(duplicateRemoverSet.begin(), duplicateRemoverSet.end());
//    duplicateRemoverSet.clear();
//
//    vector<unsigned long long int> combined(forwardFound);
//    combined.insert(combined.end(), reverseFound.begin(), reverseFound.end());
//    sort(combined.begin(), combined.end());
//    combined.erase(unique(combined.begin(), combined.end()), combined.end());

    cout << "Number of possible read locations found from forward (with duplicates): " + to_string(numLocationsForward) << endl;
//    cout << "Number of possible read locations found from forward (no duplicates): " + to_string(forwardFound.size()) << endl;
    cout << "Number of possible read locations found from reverse (with duplicates): " + to_string(numLocationsReverse) << endl;
//    cout << "Number of possible read locations found from reverse (no duplicates): " + to_string(reverseFound.size()) << endl;
    infoFile << "Number of possible read locations found from forward (with duplicates): " + to_string(numLocationsForward) << endl;
//    infoFile << "Number of possible read locations found from forward (no duplicates): " + to_string(forwardFound.size()) << endl;
    infoFile << "Number of possible read locations found from reverse (with duplicates): " + to_string(numLocationsReverse) << endl;
//    infoFile << "Number of possible read locations found from reverse (no duplicates): " + to_string(reverseFound.size()) << endl;

    cout << "Number of possible read locations from forward and reverse (with duplicates): " + to_string(numLocationsForward + numLocationsReverse) << endl << endl;
//    cout << "Number of possible read locations from forward and reverse (without duplicates): " + to_string(combined.size()) << endl;
    infoFile << "Number of possible read locations from forward and reverse (with duplicates): " + to_string(numLocationsForward + numLocationsReverse) << endl;
//    infoFile << "Number of possible read locations from forward and reverse (without duplicates): " + to_string(combined.size()) << endl;

    infoFile.close();
}

void outputPrealignmentResults() {

}