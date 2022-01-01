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
#include "inputEnums.h"
#include <limits>
#include <algorithm>
//#include "inputParser.h"
//#include <bits/stdc++.h>

using namespace std;


class Processor {
protected:
    double memorySize;
    double processorSpeed;
    Task *assignedTask;
    int assignedTaskId;
    double occupiedMemorySize;
public:
    double getOccupiedMemorySize() const;

    void setOccupiedMemorySize(double occupiedMemorySize);

public:
    bool isBusy;

    Processor() {
        this->memorySize = 0;
        this->processorSpeed = 1;
        isBusy = false;
        assignedTask = nullptr;
        assignedTaskId = -1;
    }

    explicit Processor(double memorySize) {
        this->memorySize = memorySize;
        this->processorSpeed = 1;
        isBusy = false;
        assignedTask = nullptr;
        assignedTaskId = -1;
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


    int getAssignedTaskId() const;

    void assignTaskId(unsigned int taskToBeAssigned);

    Task *getAssignedTask() const;

    void assignTask(Task *taskToBeAssigned);

    void setAssignedTaskId(int assignedTaskId);

    void setAssignedTask(Task *assignedTask);
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
    Cluster(vector<unsigned int> *groupSizes, vector<double> *memories,vector<double> *speeds) {
        if (adjacent_find( memories->begin(), memories->end(), not_equal_to<>() ) == memories->end() )
        {
            this->isMemoryHomogeneous = true;
        }
        else this->isMemoryHomogeneous = false;

        if (adjacent_find( speeds->begin(), speeds->end(), not_equal_to<>() ) == speeds->end() )
        {
            this->isProcessorHomogeneous = true;
        }
        else this->isProcessorHomogeneous = false;

        for (int i = 0; i< groupSizes->size();i++){
            unsigned int groupSize = groupSizes->at(i);
            for (unsigned int j = 0; j< groupSize; j++){
                this->processors.push_back(new Processor(memories->at(i),speeds->at(i)));
            }
        }
        initHomogeneousBandwidth(memories->size());
    }


    ~Cluster() {
        bandwidths.resize(0);

        for (unsigned long i = 0; i < processors.size(); i++) {
            delete processors.at(i);
        }
    }

    int getConfiguration() {
        //TODO: implement
        return 0;
    }

public:
    bool isHomogeneous() const {
        return isMemoryHomogeneous && isProcessorHomogeneous && isBandwidthHomogenenous;
    }

    bool isBandwidthHomogeneous() const {
        return isBandwidthHomogenenous;
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
                 << " and busy? " << (*iter)->isBusy << "assigned " << (*iter)->getAssignedTaskId() << endl;
        }
    }

    void printInfo() {
        double cumulativeMemory = 0;
        for (vector<Processor *>::iterator iter = this->processors.begin(); iter < processors.end(); iter++) {
            cumulativeMemory += (*iter)->getMemorySize();
        }
        cout << fixed << (isMemoryHomogeneous ? "Homogeneous" : "Heterogeneous") << " cluster with "
             << processors.size()
             << " processors," << " cumulative memory " << cumulativeMemory << endl;

    }

    string getPrettyClusterString() {
        string out = "Cluster:\n";
        int numProcessors = this->processors.size();
        out += "#Nodes: " + to_string(numProcessors) + ", ";
        out += "Configuration: " + to_string(this->getConfiguration()) + ", ";
        out += "MinMemory: " + to_string(this->getLastProcessorMem()) + ", ";
        out += "MaxMemory: " + to_string(this->getProcessors().at(0)->getMemorySize()) + ", ";
        out += "CumulativeMemory: " + to_string(this->getCumulativeMemory());

        /*
        out += "MinMemory: "+to_string(this->getLastProcessorMem())+", ";
        out += "CumulativeMemory: "+to_string(this->buildHomogeneousMemorySizes(maxoutd, numProcessors));*/
        return out;
    }

    int getLastProcessorMem() {
        //so far, it returns the proessor with the minimal memory size by using a linear search
        int min = std::numeric_limits<int>::max();
        for (Processor *proc: (this->processors)) {
            min = (min > proc->getMemorySize()) ? proc->getMemorySize() : min;
        }
        return min;
    }


    long getCumulativeMemory() {
        long sum = 0;
        for (Processor *proc: (this->processors)) {
            sum += proc->getMemorySize();
        }
        return sum;
    }

    string getAverageLoadAndNumberOfUsedProcessors() {
        long sumLoad = 0, sumMems = 0, numberUsed = 0;
        double percentageUsed = 0, avgLoad = 0;
        for (Processor *proc: (this->processors)) {
            if (proc->getAssignedTask() != nullptr) {
                sumLoad += proc->getOccupiedMemorySize();
                sumMems += proc->getMemorySize();
                numberUsed++;
                if (proc->getOccupiedMemorySize() == 0) {
                    cout << "task assigned and no mem occupation " << proc->getAssignedTaskId() << endl;
                }
            }
        }
        avgLoad = sumLoad * 100 / sumMems;
        percentageUsed = numberUsed * 100 / this->getProcessors().size();
        return "Load "+ to_string(avgLoad)+", Occupied "+ to_string(percentageUsed);
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

    Processor *getFirstFreeProcessorOrSmallest();

    bool hasFreeProcessor();

    void assignTasksForIds(Tree *tree);

    static void buildStatic2LevelCluster(double maxMinMem, double maxEdgesToMakespanWeights);

    static void
    buildTreeDepHomBandwidths(double CCR, unsigned int num_processors, Tree *treeobj, double &minMem, double &maxoutd,
                              schedule_traversal *&temp_schedule);

    static void buildMemHetTreeDepCluster(double CCR, unsigned int num_processors, Tree *treeobj);

    static void
    buildHomogeneousCluster(double CCR, unsigned int num_processors, Tree *treeobj, HeterogeneousAdaptationMode mode);

    static vector<double> buildNLevelMemorySizes(vector<double> memories, vector<unsigned int> processorGroupSizes);

    static void buildHomStatic2LevelCluster(double maxMinMem, double maxEdgesToMakespanWeights,
                                            HeterogeneousAdaptationMode adaptationMode);

    static void buildStatic3LevelCluster(double maxMinMem, double maxEdgesToMakespanWeights);

    static void buildHomStatic3LevelCluster(double maxMinMem, double maxEdgesToMakespanWeights,
                                            HeterogeneousAdaptationMode adaptationMode);

    string getUsageString();

    void clean();

    string getShortUsageString();

    static void
    BuildFixedClusterWithMemories(double maxEdgesToMakespanWeights, int num_processors, vector<double> &memorySizes);
};

#endif
