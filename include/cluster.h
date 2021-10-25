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

public:
    bool isHomogeneous() const {
        return isMemoryHomogeneous && isProcessorHomogeneous && isBandwidthHomogenenous;
    }

    vector<Processor *> getProcessors() {
        return this->processors;
    }

    int getNumberProcessors() {
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
        }
        this->isBandwidthHomogenenous = true;
    }

    void setMemorySizes(vector<double> &memories) {
        isMemoryHomogeneous = false;
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

    static void setFixedCluster(Cluster *cluster) {
        Cluster::fixedCluster = cluster;
    }

    static Cluster *getFixedCluster() {
        return Cluster::fixedCluster;
    }

    Processor *getFirstFreeProcessor();

    static vector<double> buildMemorySizes(double maxoutd, double minMem, unsigned int num_processors);

    static std::map<int, int> buildProcessorSpeeds(int num_processors);

    void SetBandwidth(double CCR, unsigned long tree_size, double *ewghts, double *timewghts);
};

#endif
