#include "common.h"
#include "pigeonhole.h"

string refGenome;
vector<string> reads;
map<int, vector<int>> minimizers;

int main(int argc, char *argv[]) {
    cout << "Reading the reference genome..." << endl;

    string filepath = "/home/saimanalili/multicore-rm/Pigeonhole/genome/ndna_1m.fa";

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

    int m = 100;
    int q = 8;
    int k = 2;
    int j = m / q;

    cout << "Processing minimizers..." << endl;
    string fileMinimizers = "/home/saimanalili/multicore-rm/Pigeonhole/minimizers/min_ndna_1m.txt";
    minimizers = getMinimizersFromFile(fileMinimizers);
    cout << "Done!" << endl;

    cout << "Reading the reads..." << endl;
    string readsFilename = "/home/saimanalili/multicore-rm/Pigeonhole/reads/ndna_1m_100bp_1000r_noerr.fa";
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