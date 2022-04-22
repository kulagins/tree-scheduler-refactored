//
// Created by kulagins on 22.04.22.
//
#include <string>
#include <iostream>

#ifndef TREE_SCHEDULER_REFACTORED_OUTPUTPRINTER_H
#define TREE_SCHEDULER_REFACTORED_OUTPUTPRINTER_H

using namespace std;

class OutputPrinter {
    bool verbose = true;
public:
    bool isVerbose() const {
        return this->verbose;
    }

    void setVerbose(bool verb) {
        this->verbose = verb;
    }

    void initOutput() {
        if (!verbose) {
            cout.setstate(std::ios_base::failbit);
        }
    }

    void quietPrint(string text) {
        cout.clear();
        cout << text << endl;
        this->initOutput();
    }
};

#endif
