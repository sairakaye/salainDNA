#include "bitmatrix.h"
#include "edlib.h"

//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "openmp-use-default-none"
//#include <unistd.h>



using namespace std;


int TruePos = 0;
int TrueNeg = 0;
int FalsePos = 0;
int FalseNeg = 0;

int alignmentNeeded = 0;
int notNeeded = 0;


vector<vector<int>> createMap (string P, string T, int E){
    int m = T.size();
    int n = P.size();
    char t[m];
    char p[n];
    int count;


    vector<vector<int>> NMap;
    for(int i = 0; i < 2*E + 1; i++){
        NMap.push_back(vector<int>());
    }


    int jE;
    for(int i = 0; i<m; i++){
        jE=0;
        for(int j = i-E; j<i+E+1;j++){
            if(jE <= 2*E) {
                if (T[j] != NULL && P[i] == T[j]) {
                    NMap[jE].push_back(0);
                    jE++;
                } else if (T[j] != NULL && P[i] != T[j]) {
                    NMap[jE].push_back(1);
                    jE++;
                    count++;
                }
                else{
                    jE++;
                }
            }


        }
    }

    if(count <= E){
        return vector<vector<int>>();
    }

    int E1 = E;
    int E2 = 1;

    for (int i = 0; i < E; i++){
        for (int j = 0; j < E1; j++)
        {
            NMap[i].push_back(-1);

        }
        //hammingMap[i].erase(hammingMap[i].end() - E1, hammingMap[i].end());
        E1--;
    }
    for (int i = E + 1; i < E * 2 + 1; i++){
        for (int j = 0; j < E2; j++){
            NMap[i].insert(NMap[i].begin() + j, -1);
        }
        //hammingMap[i].erase(hammingMap[i].end() - E1, hammingMap[i].end());
        E2++;

    }

    return NMap;
}


int countZeroes(vector<int> toCount){
    int ctr = 0;

    for(int i : toCount){
        if (i == 0){
            ctr++;
        }
    }
    return ctr;

}


int countOnes(vector<int> toCount, int E){
    int ctr = 0;

    if(toCount.empty()){
        return E + 1;
    }
    else {
        for (int i = 0; i < toCount.size(); i++) {
            if (toCount[i] == 1 && toCount[i] == 1) {
                ctr++;
                i++;
            } else if (toCount[i] == 1) {
                ctr++;
            }
        }
    }
    return ctr;
}
vector<int> checkDiagonals(int m, vector<vector<int>> NMap, int E, int w, int windowSize){
    vector<int> bestDiagonal = {1,1,1,1};
    vector<int> currDiagonal;

    int bestZeroCount = 0;
    int ctr = 0;

    for(int i = 0; i < E * 2 + 1; i++){
        for(int j = 0; j < windowSize; j++){
            currDiagonal.push_back(NMap[i][w+j]);
        }
        int currZeroCount = countZeroes(currDiagonal);

        if(currZeroCount >= bestZeroCount) {
            bestZeroCount = currZeroCount;
            bestDiagonal = currDiagonal;
        }
        currDiagonal.clear();
    }

    //bestDiagonalsArray[w] = bestDiagonal;
    return bestDiagonal;
}

vector<int> slidingWindow(vector<vector<int>> NMap, int m, int E) {
    int windowSize = 4;
    vector<int> finalVector;
    vector<int> mainDiagonal;
    //vector<thread> threadArray;

    mainDiagonal = NMap[E];
    finalVector = mainDiagonal;

    if(NMap.empty()){
        return vector<int>();
    }
    else {
        for (int w = 0; w < m; w++) {
            if (w + windowSize > m) {
                windowSize = m - w;
            }

            vector<int> bestDiagonalinVector = checkDiagonals(m, NMap, E, w, windowSize);

            int l = 0;
            for (int i = w; i < windowSize + w; i++) {
                if (i < m) {
                    finalVector[i] = bestDiagonalinVector[l];
                    l++;
                }
            }



/*
 *      vector<int> tempDiag;
        vector<int> seedVector;
        for (int i = w; i < windowSize + w; i++) {
            if (i < m) {
                seedVector.push_back(finalVector[i]);
                tempDiag.push_back(mainDiagonal[i]);
            }
        }


        vector<int> finalDiag;
        vector<int> tempfinalVector = finalVector;
        finalDiag = tempDiag;

        if (windowSize >= 3) {
            vector<int> finalDiag;
            finalDiag = tempDiag;

            if (countZeroes(tempDiag) < countZeroes(bestDiagonalinVector)) {
                finalDiag = bestDiagonalinVector;
            }


            if (countZeroes(seedVector) > countZeroes(finalDiag)) {
                finalDiag = seedVector;
            }

            int l = 0;
            for (int i = w; i < windowSize + w; i++) {
                if (i < m) {
                    finalVector[i] = finalDiag[l];
                    l++;
                }

            }
        }*/
        }
    }
    return finalVector;

}

vector<int> threadFunc(int E , string read, string reference){

    int m = reference.length();
    int n = read.length();


    vector<int> shouji(m);
    //vector<vector<int>> nmap;
    shouji = slidingWindow(createMap(read, reference, E),m, E);

    return shouji;
    /*if (countOnes(shouji) <= E) {
        alignmentNeeded++;
    } else {
        notNeeded++;
    }*/

}


void multiThreadedMain() {

    // string readString, referenceString;
    //double finalTime;

    //vector<string> read;
    //vector<string> reference;

    //ifstream infile("/home/aaron/Desktop/Shouji-master/Datasets/ERR240727_1_E40_30million.txt");
    //ifstream infile("/home/aaron/Desktop/10k");

    //vector<string> fileNames;
    //vector<thread> threadArray;
    cout << "Time\tE\tNeeded\tNotNeeded" << endl;

    /*
    while (infile >> readString >> referenceString) {
        read.push_back(readString);
        reference.push_back(referenceString);
    }
    */

    vector<string> reads;
    vector<string> refs;

    int refcount = refs.size();
    int readcount = reads.size();

    int E = e;

    alignmentNeeded = 0;
    notNeeded = 0;

    TruePos = 0;
    TrueNeg = 0;
    FalsePos = 0;
    FalseNeg = 0;

    auto start = std::chrono::high_resolution_clock::now();

    int i;
    #pragma omp parallel for reduction(+:notNeeded, alignmentNeeded)
    for(i = 0; i < forwardReadsMap.size(); i++) {
        auto iteratorMap = forwardReadsMap.begin();
        advance(iteratorMap, i);

        vector<unsigned long long int>& locations = (*iteratorMap).second;

        int j;
        for (j = 0; j < locations.size(); j++) {
            if (countOnes(threadFunc(E, (*iteratorMap).first, refGenome.substr(locations[j], m)), E) <= E) {
                #pragma omp critical
                filteredReadsMap[(*iteratorMap).first].push_back(locations[j]);
                alignmentNeeded++;
            } else {
                notNeeded++;
            }

        }
    }

    #pragma omp parallel for reduction(+:notNeeded, alignmentNeeded)
    for(i = 0; i < reverseReadsMap.size(); i++) {
        auto iteratorMap = reverseReadsMap.begin();
        advance(iteratorMap, i);

        vector<unsigned long long int>& locations = (*iteratorMap).second;

        int j;
        for (j = 0; j < locations.size(); j++) {
            if (countOnes(threadFunc(E, (*iteratorMap).first, refGenome.substr(locations[j], m)), E) <= E) {
                #pragma omp critical
                filteredReadsMap[(*iteratorMap).first].push_back(locations[j]);
                alignmentNeeded++;
            } else {
                notNeeded++;
            }

        }
    }

    auto end = std::chrono::high_resolution_clock::now();

    chrono::duration<double> diff = end-start;

    cout << diff.count() << "\t" <<E<<"\t"<<alignmentNeeded<<"\t"<<notNeeded << endl;
}

void checkResultswithEdlib(){


    string readString, referenceString;
    vector<string> read;
    vector<string> reference;



    int ctr = 0;

    //ifstream infile("/home/aaron/Desktop/ERR240727_1_E2_30million.txt");
    ifstream infile("/home/aaron/Desktop/Shouji-master/Datasets/ERR240727_1_E2_30million.txt");
    //ifstream infile("/home/aaron/Desktop/Shouji-master/Datasets/Test");
    //ifstream infile("F0.txt");
    while (infile >> readString >> referenceString) {
        read.push_back(readString);
        reference.push_back(referenceString);
    }

    cout << "Time\tE\tTruePos\tTrueNeg\tFalsePos\tFalseNeg";

    for (int E = 0; E < 11; E++) {
        double finalTime;
        clock_t start, end;
        TruePos = 0;
        TrueNeg = 0;
        FalsePos = 0;
        FalseNeg = 0;

        bool EdlibAccept;
        bool ShoujiAccept;

        int m = reference[0].length();

        start = clock();
        for (int i = 0; i < reference.size(); i++) {

            ctr++;
            vector<int> shouji(m);


            start = clock();

            //if (E == 4) {

            shouji = slidingWindow(createMap(read[i], reference[i], E), m, E);

            //PRINTTTT
            /* cout << "\n" << "Final bit vector : ";
             for (int i = 0; i < shouji.size(); i++) {
                 cout << shouji[i];
             }
             cout << "\n" << "NUMBER OF ONES: ";
             cout << countOnes(shouji);
             cout << "\n" << "Size: ";
             cout << shouji.size();*/

            if (countOnes(shouji, E) <= E) {
                ShoujiAccept = true;
            } else {
                ShoujiAccept = false;
            }




            const char* const edlibRead = read[i].c_str();


            EdlibAlignResult resultEdlib = edlibAlign(reference[i].c_str(), m, edlibRead, m, edlibNewAlignConfig(E, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
            edlibFreeAlignResult(resultEdlib);
            if (resultEdlib.editDistance!= -1)
                EdlibAccept = true;
            else
                EdlibAccept =false;


            if(ShoujiAccept && EdlibAccept){
                TruePos++;
            }
            else if(!ShoujiAccept && !EdlibAccept){
                TrueNeg++;
            }
            else if(ShoujiAccept && !EdlibAccept){
                FalsePos++;
            }
            else if(!ShoujiAccept && EdlibAccept){
                FalseNeg++;
            }

            //}

            //}
            end = clock();

            finalTime += double(end - start) / double(CLOCKS_PER_SEC);
        }
        cout << "\n";
        cout << finalTime << "\t" <<E<<"\t"<<TruePos<<"\t"<<TrueNeg<<"\t"<<FalsePos<<"\t"<<FalseNeg;
    }
}
/*

int main(void){
    //fixedInputMain();
    //multipleInputMain();

    multiThreadedMain();

    //checkResultswithEdlib();

    return 0;
}
*/

//#pragma clang diagnostic pop