//
// Created by saimanalili on 17/04/2020.
//

#include "output.h"

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
    string infoFileName(mode + "_info_" + mainName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
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
    cout << "Locations found by Bit Matrix: " << numFilteredReadLocations << endl;
}

void outputSAMFile() {
    ofstream SAMFile;
    SAMFile.open(SAMFileName.c_str(), ios::out);
    SAMFile << "@HD\tVN:1.4" << endl;
    SAMFile << "@SQ\tSN:" + refGenome.genomeName + "\tLN:" + to_string(refGenome.genomeData.size()) << endl;

    for (int i = 0; i < reads.size(); i++) {
        for (int j = 0; j < reads[i].forwardLocations.size(); j++) {
            EdlibAlignResult result = edlibAlign(reads[i].readData.c_str(), reads[i].readData.length(),
                    refGenome.genomeData.substr(reads[i].forwardLocations[j], m).c_str(), reads[i].readData.length(), edlibDefaultAlignConfig());
            string cigarResult(edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD));
            edlibFreeAlignResult(result);

            SAMFile << reads[i].readName << "\t" << refGenome.genomeName << "\t" << reads[i].forwardLocations[j] << "\t" << to_string(255) << "\t";
            SAMFile << cigarResult << "\t" << "*" << "0" << "\t" << "0" << "\t" << reads[i].readData << "\t" << "*" << endl;
        }

        for (int j = 0; j < reads[i].reverseLocations.size(); j++) {
            EdlibAlignResult result = edlibAlign(reads[i].readData.c_str(), m,
                                                 refGenome.genomeData.substr(reads[i].reverseLocations[j], m).c_str(), m, edlibDefaultAlignConfig());
            string cigarResult(edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD));

            SAMFile << reads[i].readName << "\t" << refGenome.genomeName << "\t" << reads[i].reverseLocations[j] << "\t";
            SAMFile << to_string(255) << cigarResult << "\t" << "*" << "0" << "\t" << "0" << reads[i].readData << "\t" << "*" << endl;
        }
    }
}