/**
 * @file output.cpp
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief It contains the implementation of the output results.
 */

#include "output.h"

/**
 * It writes a file that outputs a pair of reads from the seed selector. T
 *
 * @param genomeFileName - The filename of the reference genome.
 */
void outputPairReads(string& genomeFileName) {
    ofstream pairReadsFile;
    string pairReadsFileName(mode + "_pair_reads_" + genomeFileName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + searchMode + ".txt");
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

/**
 * It writes a file that outputs the time of each process in SalainDNA.
 *
 * @param genomeFileName - The filename of the reference genome.
 * @param indexRunTime - The runtime of the Hash-based Indexer.
 * @param ssRunTime - The runtime of Seed Selector.
 * @param bmRunTime - The runtime of the Bit Matrix.
 * @param verificationRunTime - The runtime of Edlib.
 * @param totalRunTime - The total run time of SalainDNA.
 */
void outputRunTimeResults(string& genomeFileName, double indexRunTime, double ssRunTime, double bmRunTime,
                          double verificationRunTime, double totalRunTime) {
    ofstream runTimeFile;
    string runTimeFileName(mode + "_end2end_" + genomeFileName + "_" + to_string(reads.size()) + "_" + to_string(m) + "R_" + to_string(q) + "_" + to_string(e) + "_" + searchMode + ".txt");
    runTimeFile.open(runTimeFileName.c_str(), ios::out);

    runTimeFile << "Indexing: " << to_string(indexRunTime) + " sec" << endl;
    runTimeFile << "Seed Selector: " << to_string(ssRunTime) + " sec" << endl;
    runTimeFile << "Bit Matrix: " << to_string(bmRunTime) + " sec" << endl;
    runTimeFile << "Verification: " << to_string(verificationRunTime) + " sec" << endl;
    runTimeFile << "Total time (start to bottom): " << to_string(totalRunTime) + " sec" << endl;

    runTimeFile.close();
}

/**
 * It outputs the results of the seed selector process.
 *
 * @param genomeFileName - The filename of the reference genome.
 * @param timeTaken - The run time taken by the Seed Selector.
 */
void outputSeedSelectorResults(string& genomeFileName, double timeTaken) {
    cout << "Time taken by the Seed Selector process is: " << to_string(timeTaken) << " sec" << endl << endl;

    cout << "Number of seeds checked: " << to_string(numSeeds) << endl;
    cout << "Number of reads checked: " << to_string(numReads) << endl << endl;

    cout << "Number of accepted seeds: " << to_string(numAcceptedSeeds) << endl;
    cout << "Number of accepted reads: " << to_string(numAcceptedReads) << endl << endl;

    cout << "Number of possible read locations found: " + to_string(numPossibleReadLocations) << endl << endl;
}

/**
 * It outputs the results of the pre-alignment filtering stage.
 *
 */
void outputPrealignmentResults() {
    cout << "There are " << numReads << " reads. " << numAcceptedReads << " are accepted for e = " + to_string(e) << "." << endl;
    cout << "Sequence Name: " << refGenome.genomeName << endl;
    cout << "Locations found by Seed Selector: " << numPossibleReadLocations << endl;
    cout << "Locations accepted by Bit Matrix: " << numFilteredReadLocations << endl;
}

/**
 * It outputs the results of the verification stage using Edlib.
 *
 */
void outputEdlibResults() {
    cout << "Locations accepted by Edlib: " << numVerifiedReadLocations << endl;
}

/**
 * It outputs the SAM file of the SalainDNA read mapping process.
 *
 */
void outputSAMFile() {
    ofstream SAMFile;
    SAMFile.open(SAMFileName.c_str(), ios::out);
    SAMFile << "@HD\tVN:1.4" << endl;
    SAMFile << "@SQ\tSN:" + refGenome.genomeName + "\tLN:" + to_string(refGenome.genomeData.size()) << endl;
    SAMFile << "@PG\tID:salainDNA\tPN:salainDNA\tVN:1.0" << endl;

    for (int i = 0; i < reads.size(); i++) {
        if (reads[i].forwardLocations.size() > 0) {
            string read = reads[i].readData;
            vector<unsigned long long int>& locations = reads[i].forwardLocations;

            for (int j = 0; j < locations.size(); j++) {
                EdlibAlignResult result = edlibAlign(read.c_str(), read.length(),
                                                     refGenome.genomeData.substr(locations[j], m).c_str(), read.length(),
                                                     edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                string cigarResult(edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD));
                edlibFreeAlignResult(result);

                SAMFile << reads[i].readName << "\t" << "0" << "\t" << refGenome.genomeName << "\t" << to_string(locations[j] + 1) << "\t" << to_string(255) << "\t";
                SAMFile << cigarResult << "\t" << "*" << "\t" << "0" << "\t" << "0" << "\t" << refGenome.genomeData.substr(locations[j], m) << "\t" << "*" << "\t";
                SAMFile << "NM:i:0" << "\t" << "MD:Z:" << m << endl;
            }
        }

        if (reads[i].reverseLocations.size() > 0) {
            string read = reverseComplement(reads[i].readData);
            vector<unsigned long long int>& locations = reads[i].reverseLocations;

            for (int j = 0; j < locations.size(); j++) {
                EdlibAlignResult result = edlibAlign(read.c_str(), read.length(),
                                                     refGenome.genomeData.substr(locations[j], m).c_str(), read.length(),
                                                     edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                string cigarResult(edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD));
                edlibFreeAlignResult(result);

                SAMFile << reads[i].readName << "\t" << "16" << "\t" << refGenome.genomeName << "\t" << to_string(locations[j] + 1) << "\t" << to_string(255) << "\t";
                SAMFile << cigarResult << "\t" << "*" << "\t" << "0" << "\t" << "0" << "\t" << refGenome.genomeData.substr(locations[j], m) << "\t" << "*" << "\t";
                SAMFile << "NM:i:0" << "\t" << "MD:Z:" << m << endl;
            }
        }
    }
}