//
// Created by saimanalili on 29/12/2019.
//


#ifndef MP_PIGEONHOLE_FILE_MANAGER_H
#define MP_PIGEONHOLE_FILE_MANAGER_H

#include <iostream>
#include <fstream>
#include <map>
#include <iterator>
#include <string>
#include <utility>
#include <algorithm>
#include <vector>
#include <bits/stdc++.h>

typedef struct {
    char *genome;
    int length;
} RefGenome;

typedef struct {
    char **reads;
    int size;
} ReadList;

std::string getFilename(std::string filename);
RefGenome *readGenome(std::string filename);
ReadList *readReads(std::string filename);
std::map<int, std::vector<int>> getMinimizersFromFile(std::string filename);

#endif //MP_PIGEONHOLE_FILE_MANAGER_H
