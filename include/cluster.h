#ifndef cluster_h
#define cluster_h

#include <iostream>
#include <list>
#include <ostream>
#include <stdio.h>
#include <assert.h>
#include <forward_list>
#include <map>
#include <vector>
#include "tree.fwd.h"

using namespace std;

class Processor {
protected:
    double memorySize;
    double processorSpeed;
    Task *assignedTask;

public:
    bool isBusy;

    Processor() {
        this->memorySize = 0;
        this->processorSpeed = 1;
        isBusy = false;
        assignedTask = nullptr;
    }

    explicit Processor(double memorySize) {
        this->memorySize = memorySize;
        this->processorSpeed = 1;
        isBusy = false;
        assignedTask = nullptr;
    }

    Processor(double memorySize, double processorSpeed) {
        this->memorySize = memorySize;
        this->processorSpeed = processorSpeed;
        isBusy = false;
        assignedTask = nullptr;
    }

    double getMemorySize() const {
        return memorySize;
    }

    void setMemorySize(double memory) {
        this->memorySize = memory;
    }

    double getProcessorSpeed() const {
        return processorSpeed;
    }

    void assignTask(Task *taskToBeAssigned) {
        this->assignedTask = taskToBeAssigned;
        this->isBusy = true;
    }
};

class Cluster {
protected:
    bool isMemoryHomogeneous;
    bool isProcessorHomogeneous;
    bool isBandwidthHomogenenous;

    vector<Processor *> processors;
    vector<vector<double>> bandwidths;
    static Cluster *fixedCluster;

public:
    Cluster() {
        this->isMemoryHomogeneous = this->isProcessorHomogeneous = this->isBandwidthHomogenenous = true;

        processors.resize(3, new Processor());
        bandwidths.resize(3);
        for (unsigned long i = 0; i < bandwidths.size(); i++) {
            bandwidths.at(i).resize(3, 1);
        }
    }

    Cluster(unsigned int clusterSize, bool isMemoryHomogeneous) {
        this->isMemoryHomogeneous = isMemoryHomogeneous;
        this->isProcessorHomogeneous = this->isBandwidthHomogenenous = true;

        processors.resize(clusterSize);
        for (unsigned long i = 0; i < clusterSize; i++) {
            processors.at(i) = new Processor();
        }
        initHomogeneousBandwidth(clusterSize);
    }

    Cluster(vector<double> memories) {
        isProcessorHomogeneous = isBandwidthHomogenenous = true;
        setMemorySizes(memories);
        initHomogeneousBandwidth(memories.size());
    }


    ~Cluster() {
        bandwidths.resize(0);

        for (unsigned long i = 0; i < processors.size(); i++) {
            delete processors.at(i);
        }
    }

    int getConfiguration(){
        //TODO: implement
        return 0;
    }
public:
    bool isHomogeneous() const {
        return isMemoryHomogeneous && isProcessorHomogeneous && isBandwidthHomogenenous;
    }

    void setHomogeneity(bool homogeneity) {
        isMemoryHomogeneous = homogeneity;
    }

    vector<Processor *> getProcessors() {
        return this->processors;
    }

    unsigned int getNumberProcessors() {
        return this->processors.size();
    }

    double getBandwidthBetween(int firstProcessor, int secondProcessor) {
        return this->bandwidths.at(firstProcessor).at(secondProcessor);
    }

    double getHomogeneousBandwidth() {
        return this->bandwidths.at(0).at(1);
    }

    void initHomogeneousBandwidth(int bandwidthsNumber, double bandwidth = 1) {
        bandwidths.resize(bandwidthsNumber);
        setHomogeneousBandwidth(bandwidth);
    }

    void setHomogeneousBandwidth(double bandwidth) {
        for (unsigned long i = 0; i < bandwidths.size(); i++) {
            //TODO init only upper half
            bandwidths.at(i).resize(bandwidths.size(), bandwidth);
            for (unsigned long j = 0; j < bandwidths.size(); j++) {
                bandwidths.at(i).at(j) = bandwidth;
            }
        }
        this->isBandwidthHomogenenous = true;
    }

    void setMemorySizes(vector<double> &memories) {
        //  isMemoryHomogeneous = false;
        if (processors.size() != memories.size()) {
            processors.resize(memories.size());
            for (unsigned long i = 0; i < memories.size(); i++) {
                processors.at(i) = new Processor(memories.at(i));
            }
        } else {
            for (unsigned long i = 0; i < memories.size(); i++) {
                processors.at(i)->setMemorySize(memories.at(i));
            }
        }

    }

    void printProcessors() {
        for (vector<Processor *>::iterator iter = this->processors.begin(); iter < processors.end(); iter++) {
            cout << "Processor with memory " << (*iter)->getMemorySize() << ", speed " << (*iter)->getProcessorSpeed()
                 << " and busy? " << (*iter)->isBusy << endl;
        }
    }

    string getPrettyClusterString(){
        string out = "Cluster:\n";
        int numProcessors = this->processors.size();
        out += "#Nodes: "+to_string(numProcessors)+", ";
        out += "Configuration: "+to_string(this->getConfiguration())+", ";
        out += "MinMemory: "+to_string(this->getLastProcessor())+", ";
        out += "CumulativeMemory: "+to_string(this->getCumulativeMemory());
        /*
        out += "MinMemory: "+to_string(this->getLastProcessor())+", ";
        out += "CumulativeMemory: "+to_string(this->buildHomogeneousMemorySizes(maxoutd, numProcessors));*/
        return out;
    }

    int getLastProcessor(){
        //so far, it returns the proessor with the minimal memory size by using a linear search
        int min = INT_MAX;
        for (Processor *proc: (this->processors)){
            min = (min>proc->getMemorySize())? proc->getMemorySize():min;
        }
        return min;
    }

    int getCumulativeMemory(){
        int sum = 0;
        for (Processor *proc: (this->processors)){
            sum += proc->getMemorySize();
        }
        return sum;
    }

    static void setFixedCluster(Cluster *cluster) {
        Cluster::fixedCluster = cluster;
    }

    static Cluster *getFixedCluster() {
        return Cluster::fixedCluster;
    }

    Processor *getFirstFreeProcessor();

    static vector<double> build3LevelMemorySizes(double maxoutd, double minMem, unsigned int num_processors);

    static vector<double> buildHomogeneousMemorySizes(double memSize, unsigned int num_processors);

    static std::map<int, int> buildProcessorSpeeds(int num_processors);

    void SetBandwidth(double CCR, Tree *treeobj);
};

#endif
