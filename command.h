/**
 * @file command.h
 * @author A. Fajardo, S. Manalili, C. Mercado, R. Zapanta
 * @brief Header file for processing the Commands, containing function declaration to be used by the program.
 */

#ifndef SALAINDNA_COMMAND_H
#define SALAINDNA_COMMAND_H

#include "common.h"

void processingArguments(int argc, char *argv[], string &genomeFilePath, string &readsFilePath, string &indexFilePath, string &genomeFileName);

#endif
