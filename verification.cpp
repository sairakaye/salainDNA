//
// Created by saimanalili on 21/09/2020.
//

#include "verification.h"

void verifyWithEdlib() {
    //cout << "Edlib:\nTime\tE\tAccepted\tRejected" << endl;

    auto start = std::chrono::high_resolution_clock::now();
    int i;
    for (i = 0; i < reads.size(); i++) {
        string read = reads[i].readData;
        //vector<unsigned long long int> tempAcceptedLocations;

        if (reads[i].forwardLocations.size() > 0) {
            //vector<unsigned long long int> &locations = reads[i].forwardLocations;

            int j;
            for (j = 0; j < reads[i].forwardLocations.size(); j++) {
                if (refGenome.genomeData.substr(reads[i].forwardLocations[j], m).size() == m) {
                    EdlibAlignResult resultEdlib;
                    string refGenomeRead = refGenome.genomeData.substr(reads[i].forwardLocations[j], m);
                    const char* const pRef = refGenomeRead.c_str();
                    const char* const pRead = read.c_str();
                    resultEdlib = edlibAlign(pRef, m, pRead, m,
                                             edlibNewAlignConfig(e, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                    edlibFreeAlignResult(resultEdlib);
                    if (resultEdlib.editDistance != -1) {
                        //tempAcceptedLocations.push_back(locations[j]);
                        truePos++;
//                        EdlibAccept = true;
                    } else {
//                        EdlibAccept = false;
                        trueNeg++;
                        reads[i].forwardLocations.erase(reads[i].forwardLocations.begin() + j);
                        j--;

                    }
                } else {
                    trueNeg++;
                    reads[i].forwardLocations.erase(reads[i].forwardLocations.begin() + j);
                    j--;
                }
//                auto end = std::chrono::high_resolution_clock::now();
            }
            //reads[i].forwardLocations = vector<unsigned long long int>(tempAcceptedLocations);
        }

        if (reads[i].reverseLocations.size() > 0) {
            //vector<unsigned long long int> &locations = reads[i].reverseLocations;
            string reverseRead = reverseComplement(read);

            int j;
            for (j = 0; j < reads[i].reverseLocations.size(); j++) {
                if (refGenome.genomeData.substr(reads[i].reverseLocations[j], m).size() == m) {
                    EdlibAlignResult resultEdlib;
                    string refGenomeRead = refGenome.genomeData.substr(reads[i].reverseLocations[j], m);
                    const char* const pRef = refGenomeRead.c_str();
                    const char* const pRead = reverseRead.c_str();
                    resultEdlib = edlibAlign(pRef, m, pRead, m,
                                             edlibNewAlignConfig(e, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                    edlibFreeAlignResult(resultEdlib);
                    if (resultEdlib.editDistance!= -1) {
                        //tempAcceptedLocations.push_back(locations[j]);
                        truePos++;
//                        EdlibAccept = true;
                    } else {
//                        EdlibAccept = false;
                        trueNeg++;
                        reads[i].reverseLocations.erase(reads[i].reverseLocations.begin() + j);
                        j--;
                    }
                } else {
                    trueNeg++;
                    reads[i].reverseLocations.erase(reads[i].reverseLocations.begin() + j);
                    j--;
                }
//                auto end = std::chrono::high_resolution_clock::now();
            }
            //reads[i].reverseLocations = vector<unsigned long long int>(tempAcceptedLocations);
        }
    }
    numVerifiedReadLocations = truePos;

    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;
    //cout << diff.count() << "\t" << e << "\t" << truePos << "\t" << trueNeg << endl << endl;

    verificationRunTime = diff.count();
    cout << "Time taken by Edlib is: " << to_string(verificationRunTime) << " sec" << endl;
    cout << "Number of accepted locations: " << to_string(truePos) << endl;
    cout << "Number of rejected locations: " << to_string(trueNeg) << endl << endl;
}
