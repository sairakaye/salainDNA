//
// Created by saimanalili on 08/01/2020.
//

#include "minimizers.h"
using namespace std;

vector<pair<string, int> > alphabetRef = { {"A", 0}, {"C",1}, {"G",2}, {"T", 3} };

int extractRanking(string kmer)
{
    string binary;
    int rankValue;

    for (int i = 0; i<kmer.length(); i++){
        for (int j = 0; j< alphabetRef.size(); j++){
            if (kmer.at(i) + string() == alphabetRef.at(j).first){
                rankValue = alphabetRef.at(j).second;
                binary.append(bitset<2>(rankValue).to_string());
            }
        }
    }

    int decimal = stoi (binary, nullptr, 2);
    return decimal;
}