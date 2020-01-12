#include "common.h"
#include "pigeonhole.h"

string refGenome;
vector<string> reads;
map<int, vector<int>> minimizers;

int main(int argc, char *argv[]) {
    cout << "Reading the reference genome..." << endl;

    string filepath = "sim_255.fa";

    filepath.erase(remove(filepath.begin(), filepath.end(), '\"'), filepath.end());
    filepath.erase(remove(filepath.begin(), filepath.end(), '\''), filepath.end());

    refGenome = readGenome(filepath);

    cout << "Done!" << endl;

    /**
     * m - length of the read.
     * q - length of the seed (k-mer/q-gram)
     * k - errors
     * j - number of partitions in the read.
     */

    int m = 120;
    int q = 8;
    int k = 2;
    int j = m / q;

    cout << "Processing minimizers..." << endl;
    string fileMinimizers = "min_255.txt";
    minimizers = getMinimizersFromFile(fileMinimizers);
    cout << "Done!" << endl;

    cout << "Reading the reads..." << endl;
    string readsFilename = "sim_255_120bp_10r_perf.fa";
    reads = readReads(readsFilename);
    cout << "Done!" << endl;

    long start = clock();
    cout << "Doing pigeonhole process..." << endl;
    filterReads(readsFilename, q, m, j, k);
    //parallelizeFilterReads(readsFilename, q, m, j, k);
    long end = clock();
    cout << "Done!" << endl;

    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by the pigeonhole process is : " << fixed
         << time_taken << " sec" << endl;

    return 0;
}