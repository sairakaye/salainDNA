//
// Created by saimanalili on 17/04/2020.
//

#include "output.h"

void outputPairReads(string& mainName) {
    ofstream pairReadsFile;
    string pairReadsFileName(mode + "_pair_reads_" + mainName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
    pairReadsFile.open(pairReadsFileName.c_str(), ios::out);

    for (int i = 0; i < reads.size(); i++) {
        string mainRead = reads[i].readData;

        vector<unsigned long long int>& forwardLocations = reads[i].forwardLocations;
        for (int j = 0; j < forwardLocations.size(); j++) {
            pairReadsFile << mainRead << "\t" << refGenome.genomeData.substr(forwardLocations[j], m) << endl;
        }

        vector<unsigned long long int>& reverseLocations = reads[i].reverseLocations;
        for (int j = 0; j < reverseLocations.size(); j++) {
            pairReadsFile << mainRead << "\t" << refGenome.genomeData.substr(reverseLocations[j], m) << endl;
        }
    }

    pairReadsFile.close();
}

void outputRunTimeResults(string& mainName, double indexTimeTaken, double ssTimeTaken, double bmTimeTaken) {
    ofstream runTimeFile;
    string runTimeFileName(mode + "_end2end_" + mainName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + to_string(e) + searchMode + ".txt");
    runTimeFile.open(runTimeFileName.c_str(), ios::out);

    runTimeFile << "Indexing: " << to_string(indexTimeTaken) + " sec" << endl;
    runTimeFile << "Seed Selector: " << to_string(ssTimeTaken) + " sec" << endl;
    runTimeFile << "Bit Matrix: " << to_string(bmTimeTaken) + " sec" << endl;

    runTimeFile.close();
}

void outputSeedSelectorResults(string& mainName, double timeTaken) {
    cout << "Time taken by the pigeonhole process is: " << to_string(timeTaken) << " sec" << endl << endl;

    cout << "Number of seeds checked: " << to_string(numSeeds) << endl;
    cout << "Number of reads checked: " << to_string(numReads) << endl << endl;

    cout << "Number of accepted seeds: " << to_string(numAcceptedSeeds) << endl;
    cout << "Number of accepted reads: " << to_string(numAcceptedReads) << endl << endl;

    cout << "Number of possible read locations found: " + to_string(numPossibleReadLocations) << endl << endl;
}

void outputFileSeedSelectorResults(string& mainName, double timeTaken) {
    ofstream infoFile;
    string infoFileName(mode + "_info_" + mainName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + to_string(e) + searchMode + ".txt");
    infoFile.open(infoFileName.c_str(), ios::out);

    infoFile << "Time taken by the pigeonhole process is : " << to_string(timeTaken) << " sec" << endl << endl;

    infoFile << "Number of seeds checked: " << to_string(numSeeds) << endl;
    infoFile << "Number of reads checked: " << to_string(numReads) << endl << endl;

    infoFile << "Number of accepted seeds: " << to_string(numAcceptedSeeds) << endl;
    infoFile << "Number of accepted reads: " << to_string(numAcceptedReads) << endl << endl;

    infoFile << "Number of possible read locations found from forward: " + to_string(numPossibleReadLocations) << endl;

    infoFile.close();
}

void outputPrealignmentResults() {
    cout << "There are " << numReads << " reads. " << numAcceptedReads << " are accepted for e = " + to_string(e) << "." << endl;
    cout << "Sequence Name: " << refGenome.genomeName << endl;
    cout << "Locations found by Seed Selector: " << numPossibleReadLocations << endl;
    cout << "Locations accepted by Bit Matrix: " << numFilteredReadLocations << endl;
}

void outputSAMFile() {
    ofstream SAMFile;
    SAMFile.open(SAMFileName.c_str(), ios::out);
    SAMFile << "@HD\tVN:1.4" << endl;
    SAMFile << "@SQ\tSN:" + refGenome.genomeName + "\tLN:" + to_string(refGenome.genomeData.size()) << endl;
    SAMFile << "@PG\tID:multicore-rm\tPN:multicore-rm\tVN:1.0" << endl;

    for (int i = 0; i < reads.size(); i++) {
        if (reads[i].forwardLocations.size() > 0) {
            string read(reads[i].readData);
            vector<unsigned long long int>& locations = reads[i].forwardLocations;

            for (int j = 0; j < locations.size(); j++) {
                EdlibAlignResult result = edlibAlign(read.c_str(), read.length(),
                                                     refGenome.genomeData.substr(locations[j], m).c_str(), read.length(),
                                                     edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                string cigarResult(edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD));
                edlibFreeAlignResult(result);

                SAMFile << reads[i].readName << "\t" << "0" << "\t" << refGenome.genomeName << "\t" << to_string(locations[j] + 1) << "\t" << to_string(255) << "\t";
                SAMFile << cigarResult << "\t" << "*" << "\t" << "0" << "\t" << "0" << "\t" << refGenome.genomeData.substr(locations[j], m) << "\t" << "*" << endl;
            }
        } else if (reads[i].reverseLocations.size() > 0) {
            string read(reads[i].readData);
            vector<unsigned long long int>& locations = reads[i].reverseLocations;

            for (int j = 0; j < locations.size(); j++) {
                EdlibAlignResult result = edlibAlign(read.c_str(), read.length(),
                                                     refGenome.genomeData.substr(locations[j], m).c_str(), read.length(),
                                                     edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                string cigarResult(edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD));
                edlibFreeAlignResult(result);

                SAMFile << reads[i].readName << "\t" << "16" << "\t" << refGenome.genomeName << "\t" << to_string(locations[j] + 1) << "\t" << to_string(255) << "\t";
                SAMFile << cigarResult << "\t" << "*" << "\t" << "0" << "\t" << "0" << "\t" << refGenome.genomeData.substr(locations[j], m) << "\t" << "*" << endl;
            }
        }
    }
}